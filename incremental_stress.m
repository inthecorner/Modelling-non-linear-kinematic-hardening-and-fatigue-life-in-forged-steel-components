clear all
close all

% Path to output file destination
outpath = 'OUTDATA/';

% Name of output file (no suffix)
outname = 'incremental_stress';

% Calculate principal directions every increment (1) (recommended for 
% non-proportional loading) or only from the first increment (0)? 
update_principal = 0;

% =========================
%% (INPUT) DEFINE MATERIAL
% =========================
% Basic information
E   = 200e3;	% Youngs modulus [MPa]
nu  = 0.3;		% Poissons ratio [-]

% Hardening model
% Options:
%		- isotropic
%		- prager   (kinematic)
%		- ziegler  (kinematic)
%		- chaboche (combined)
hardtype = 'prager';

% ======== MATERIAL DATA FOR ISOTROPIC/LINEAR KINEMATIC MODELS =========
% Hardening data for isotropic and linear kinematic models (Abaqus format)
%								 Yield stress [MPa]		  Plastic strain [-]
matparams.tab = [350										0
								 390										0.001];

% ================ MATERIAL DATA FOR CHABOCHE MODEL ====================
% Linear term parameter [MPa]
matparams.chab.c	 = [50000		% Backstress 1
											100000];% Backstress 2

% Recall term parameter [-]
matparams.chab.gamma = [200		% Backstress 1
												2000];% Backstress 2
                      
% Isotropic hardening parameters
matparams.chab.b  = 0.2;		% Rate parameter [-]
matparams.chab.Q  = -100;		% Scale parameter [MPa]
matparams.chab.sy = 350;    % Initial yield [MPa]

% Isotropic hardening model
% Options:
%		1 - Voce (traditional)
%		2 - Hollomon (not recommended)
%		3 - Modified Hollomon (recommended)
matparams.chab.hardtype = 3;

% =========================
%% (INPUT) DEFINE LOADING
% =========================
% Number of increments per step
res = 1000;

% Define load points by target stress. Initial row must be all zero.
%						|---- Normal -----|  |------- Shear------|
%						s_11	s_22		s_33		s_12		s_13		s_23	  resolution
loadpoints=[0			0				0				0				0				0				0 	  % Initial
					  660		0				0				0				0				0				res		% Step 1
					 -450		0				0				0				0				0				res		% Step 2
					  660		0				0				0				0				0				res		% Step 3
					 -450		0				0				0				0				0				res		% Step 4
					  660		0				0				0				0				0				res];	% Step 5

% =========================
%% LOOP SETUP
% =========================

% Ensure proper string format
hardtype = lower(hardtype);

% Define loop variables
step		= 1;				% Inital step
n_steps	= size(loadpoints,1)-1; % Number of analysis steps
nits		= 0;				% Number of iterations
skip		= 0;				% Scalar for stress increment correction
dif			= zeros(3); % Scalar for stress increment correction

% Data processing on material model
if strcmp(hardtype,'chaboche')
	% Get number of backstress components
	n_alpha = length(matparams.chab.c);
else
	% Calculate plastic tangent modulus
	for i=1:size(matparams.tab,1)-1
		dy=matparams.tab(i,1) - matparams.tab(i+1,1);
		dx=matparams.tab(i,2) - matparams.tab(i+1,2);
		matparams.tab(i,3) = dy/dx;
	end
	matparams.tab(end,3) = matparams.tab(end-1,3);
	% No backstress components, only one total
	n_alpha = 0;
end

% Calculate load increments
loadvec = zeros(3,3,n_steps-1);
incvec  = zeros(3,3,n_steps-1);
for i=2:n_steps+1
	loadvec(:,:,i-1) = [loadpoints(i,[1,4,5]);
											loadpoints(i,[4,2,6]);
											loadpoints(i,[5,6,3])];
	incvec(:,:,i-1)  = [diff(loadpoints(i-1:i,[1,4,5]));
											diff(loadpoints(i-1:i,[4,2,6]));
											diff(loadpoints(i-1:i,[5,6,3]))]/loadpoints(i,7);
end

% Initial yield stress
if logical(sum(strcmp(hardtype,{'isotropic','prager','ziegler'})))
    s_y	  = matparams.tab(1,1);
elseif strcmp(hardtype,'chaboche')
    s_y	  = matparams.chab.sy;
end

% Initialise variables
s_max				= s_y;			      % Stress history internal variable
p						= 0;						  % Accumulated equivalent plastic strain
p_dot				= 0;						  % Equivalent plastic strain increment 
R						= 0;						  % Expansion of yield surface
R_dot				= 0;						  % Expansion increment of yield surface
s						= zeros(3);			  % Effective stress tensor
s_ref				= zeros(3);			  % Total stress tensor
s_dev				= zeros(3);			  % Effective deviatoric stress
s_e					= 0;						  % Equivalent effective stress
s_e_dot			= 0;						  % Equivalent effective stress increment
alpha				= zeros(3,3,n_alpha+1);	% Backstress
alpha_dot		= zeros(3,3,n_alpha+1);	% Backstress increment
e						= zeros(3);				% Total strain
e_e					= zeros(3);				% Elastic strain component
e_p					= zeros(3);				% Plastic strain component
e_p_dot			= zeros(3);				% Plastic strain increment

% Guess number of analysis points
guess = sum(loadpoints(:,7));

% Initialise history arrays for plotting
Z    = zeros(guess,1);
hist = struct;
% Absolute states
hist.s11 = Z; hist.e11 = Z; hist.ee11 = Z; hist.a11 = Z;
hist.s22 = Z; hist.e22 = Z; hist.ee22 = Z; hist.a22 = Z;
hist.s33 = Z; hist.e33 = Z; hist.ee33 = Z; hist.a33 = Z;
hist.s12 = Z; hist.e12 = Z; hist.ee12 = Z; hist.a12 = Z;
hist.s13 = Z; hist.e13 = Z; hist.ee13 = Z; hist.a13 = Z;
hist.s23 = Z; hist.e23 = Z; hist.ee23 = Z; hist.a23 = Z;
% Increment states
hist.sdot11 = Z; hist.edot11 = Z; hist.eedot11 = Z; hist.adot11 = Z;
hist.sdot22 = Z; hist.edot22 = Z; hist.eedot22 = Z; hist.adot22 = Z;
hist.sdot33 = Z; hist.edot33 = Z; hist.eedot33 = Z; hist.adot33 = Z;
hist.sdot12 = Z; hist.edot12 = Z; hist.eedot12 = Z; hist.adot12 = Z;
hist.sdot13 = Z; hist.edot13 = Z; hist.eedot13 = Z; hist.adot13 = Z;
hist.sdot23 = Z; hist.edot23 = Z; hist.eedot23 = Z; hist.adot23 = Z;
% Principal
hist.s1 = Z; hist.edot1 = Z; hist.a1 = Z;
hist.s2 = Z; hist.edot2 = Z; hist.a2 = Z;
hist.s3 = Z; hist.edot3 = Z; hist.a3 = Z;
% Misc
hist.smax = Z; hist.step = Z; hist.R = Z; hist.Rdot = Z; 
hist.p = Z; hist.pdot = Z; hist.se = Z; hist.ae = Z;

% =========================
%% MAIN LOOP
% =========================
while step <= n_steps
	% Determine the stress increment
	if skip~=1
		nits = nits+1;
		if all(dif == 0)
			s_dot = incvec(:,:,step);
		else
			s_dot = dif;
			dif = zeros(3);
		end
	end
	skip=0;
	
	% Hydrostatic tensors
	s_hyd			= eye(3)*trace(s)/3;
	s_dot_hyd = eye(3)*trace(s_dot)/3;

	% Deviatoric tensors
	s_dev			= s     - s_hyd;
	s_dot_dev = s_dot - s_dot_hyd;

	% Equivalent scalars (Mises)
	s_e	= sqrt(3/2*trace(s_dev*s_dev'));
	if s_e == 0
		s_e_dot = sqrt(3/2*trace(s_dot_dev*s_dot_dev'));
	else
		s_e_dot = 3/2*trace(s_dev*s_dot')/s_e;
	end
	
	% Make sure consistency condition is met (land on yield surface)
	if s_max-s_e>1e-9 && (s_e+s_e_dot)-s_max>1e-9
		% Calculate stress increment that lands us on surface
		s_dot_mod = s_dot*(s_max-s_e)/s_e_dot;
		% Save the excess stress in dif variable to be applied next iteration
		dif = s_dot-s_dot_mod;
		% Correct the stress increment and jump back to top of script
		s_dot = s_dot_mod;
		skip  = 1;
		continue
	elseif s_e-s_max > 1e-9 & s_e_dot>0
		s_max=s_e;
	end
	
	% Plastic switch
	if abs(s_e-s_max)<1e-9 && (s_e_dot > 0)
		% Yield surface normal	
		N = 3/2*s_dev/s_e;
		switch hardtype
		% ========= ISOTROPIC HARDENING ============
		case 'isotropic'
			% Calculate plastic tanget modulus
			h = interp1(matparams.tab(:,2),matparams.tab(:,3),p, ...
				  "previous","extrap");
			% New max stress
			R_dot = s_e_dot;

		% ========= PRAGER KINEMATIC HARDENING ============
		case 'prager'
			% Calculate plastic tanget modulus
			h = interp1(matparams.tab(:,2),matparams.tab(:,3),p, ...
				  "previous","extrap");
			% Backstress increment
			alpha_dot = s_dev*trace(N*s_dot')/s_e;

		% ========= ZIEGLER KINEMATIC HARDENING ============
		case 'ziegler'
			% Calculate plastic tanget modulus
			h = interp1(matparams.tab(:,2),matparams.tab(:,3),p, ...
				  "previous","extrap");
			% Backstress increment
			alpha_dot = s*trace(N*s_dot')/s_e;

		% ========= CHABOCHE COMBINED HARDENING ============
		case 'chaboche'
			% Model constants/parameters
			c = matparams.chab.c; gamma	= matparams.chab.gamma;
			b	= matparams.chab.b; Q	= matparams.chab.Q;
			C1 = sum(c); C2 = sum(alpha(:,:,2:end).*reshape(gamma,1,1,[]),3);
			
			% For isotropic hardening, get R_dot/p_dot value (C3)
      if matparams.chab.hardtype == 1			% Voce
				C3	= b*(Q-R);
			elseif matparams.chab.hardtype == 2	% Hollomon
				C3	= R*b/((R/Q)^(1/b));
			elseif matparams.chab.hardtype == 3	% Modified Hollomon
				v		= 10*(abs(Q)*b/sum(c))^(1/(1-b));
				C3	= (R+Q*v^b)*b/(((R+Q*v^b)/Q)^(1/b));
      end

			% Hardening modulus
			h = C3+C1-trace(C2*N');
			% Eq. plastic strain increment (lambda_dot = p_dot)
			lambda_dot = trace(N*s_dot')/h;
			% Hardening component increment
      R_dot = C3*lambda_dot;

			% Backstress increments
			alpha_dot(:,:,1) = lambda_dot*(2/3*N*C1-C2);
			for i=1:n_alpha
				alpha_dot(:,:,i+1) = lambda_dot*(2/3*N*c(i)-...
													   gamma(i)*alpha(:,:,i+1)); 
			end
		end
		% Calculate plastic strain increment
		e_p_dot = 1/h*trace(N*s_dot')*N;
	else
		% No plastic development, only elastic
		e_p_dot			= zeros(3);
		alpha_dot		= alpha*0;
		R_dot				= 0;
    N						= nan(3);
		h						= nan;
	end
	
	% Calculate equivalent plastic strain increment
	p_dot  = sqrt(2/3*trace(e_p_dot*e_p_dot'));
  
	% Calculate the pseudo-time rate of the yield surface normal
	alpha_dev_dot = alpha_dot(:,:,1)-1/3*eye(3)*trace(alpha_dot(:,:,1));
  s_dev_dot_eff = s_dot-1/3*eye(3)*trace(s_dot)-alpha_dev_dot;
  N_dot					= 1/s_e*(3/2*s_dev_dot_eff-N*trace(N*s_dev_dot_eff'));
  
	% Calculate the hardening, geometrical and total tangent moduli
	h_n		= h;
  h_g		= trace(alpha(:,:,1)*N_dot')/p_dot;
  h_tot = h_n + h_g;
	
	% Elastic strain increment
	e_e_dot = stress2estrain(E,nu,s_dot);

	% Advance variable states
	s_ref = s_ref + s_dot;									% Total stress
	alpha = alpha + alpha_dot;							% Backstress
	s     = s_ref - alpha(:,:,1);						% Effective stress
	s_dev = s-eye(3)*1/3*trace(s);					% Deviatoric effective stress
	s_e   = sqrt(3/2*trace(s_dev*s_dev'));	% Equivalent effective stress
	R     = R   + R_dot;										% Isotropic hardening variable
	s_max = s_y + R;												% Yield surface size
	p			= p   + p_dot;										% Accumulated plastic strain
	e_p	  = e_p + e_p_dot;									% Plastic strain
	e_e	  = e_e + e_e_dot;									% Elastic strain
	e     = e_e	+ e_p;											% Total strain
	a_e   = trace(alpha(:,:,1)*N');					% Equivalent backstress
  sigma = s_e + a_e;											% Equivalent total stress

	% Calculate principal stresses and strains
	if nits==1 | update_principal == 1
		[vecs,vals] = eig(s);
		[~,order]=sort(diag(vals),'descend');
		T=vecs(:,order);
	end
	principal					= diag(inv(T)*s_ref*T);
	alpha_principal		= diag(inv(T)*alpha(:,:,1)*T);
	pstrain_principal	= diag(inv(T)*e_p_dot*T);
	
	% Add tensors and scalars to history arrays
	hist = stackMat(hist, s_ref, 's', nits);
	hist = stackMat(hist, e_p,'e', nits);
	hist = stackMat(hist, e_e,'ee', nits);
	hist = stackMat(hist, e_e+e_p,'etot', nits);
	hist = stackMat(hist, alpha(:,:,1), 'a',nits);
	hist = stackMat(hist, s_dot, 'sdot', nits);
	hist = stackMat(hist, e_p_dot, 'edot', nits);
	hist = stackMat(hist, e_e_dot, 'eedot', nits);
	hist = stackMat(hist, alpha_dot(:,:,1), 'adot', nits);
	hist = stackVec(hist, principal, 'sp', nits);
	hist = stackVec(hist, alpha_principal, 'ap', nits);
	hist = stackVec(hist, pstrain_principal, 'edotp', nits);
	hist = stackScal(hist,[s_max,step,R,R_dot,p,p_dot,s_e,a_e,sigma,h_n,...
		               h_g,h_tot],{'smax', 'step', 'R', 'Rdot', 'p','pdot',...
									 'se','ae','sigma','hn','hg','htot'},nits);

	% If end load is reached, switch step
	if all( abs(loadvec(:,:,step) - s_ref) < 1e-9 ,'all')
		step=step+1;
	end
end

% Save analysis name
hist.name = outname;

% Save history to .mat file
save([outpath,outname,'.mat'],'-struct','hist');

% =========================
%% HELPER FUNCTIONS 
% =========================
function hist = stackMat(hist,matIn,name,nits)
	hist.([name,'11'])(nits) = matIn(1,1);
	hist.([name,'22'])(nits) = matIn(2,2);
	hist.([name,'33'])(nits) = matIn(3,3);
	hist.([name,'12'])(nits) = matIn(1,2);
	hist.([name,'13'])(nits) = matIn(1,3);
	hist.([name,'23'])(nits) = matIn(2,3);
end
function hist = stackVec(hist,vecIn,name,nits)
	hist.([name,'1'])(nits) = vecIn(1);
	hist.([name,'2'])(nits) = vecIn(2);
	hist.([name,'3'])(nits) = vecIn(3);
end
function hist = stackScal(hist,scalIn,names,nits)
	for i=1:length(names) 
		hist.(names{i})(nits) = scalIn(i);
	end
end
function e_e_dot = stress2estrain(E,nu,s_dot)
	d = eye(3);
	e_e_dot = zeros (3,3);
	for i1=1:3
		for j1=1:3
			sum0=0;
			for k1=1:3
				for l1=1:3
					sum0=sum0+1/E*((1+nu)/2*(d(i1,k1)*d(j1,l1)+d(i1,l1)*d(j1,k1))-...
							 nu*d(i1,j1)*d(k1,l1))*s_dot(k1,l1);
				end
			end
			e_e_dot(i1,j1) = sum0;
		end
	end
end


