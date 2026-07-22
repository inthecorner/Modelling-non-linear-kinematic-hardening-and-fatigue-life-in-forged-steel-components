function hist = total_form(s_y,E,n_alpha,c,gamma,Q,b,load,hardtype,outname)
% ----- INPUT PARAMETERS --------
%	s_y				-	Initial yield stress [MPa]
%	E					- Young's modulus [MPa]
%	n_alpha		-	Number of back-stresses
%	c					-	1D array with c values [MPa], length must match n_alpha
%	gamma			- 1D array with gamma values, length must match n_alpha
%	Q					- Isotropic hardening scaling parameter [MPa]
%	b					-	Isotropic hardening rate parameter
%	load			-	Cell array with each entry containing strain sample points
% hardtype	- Type of isotropic hardening
% outname		- Name to save the analysis under

% Ensure integer format 
n_alpha  = round(n_alpha);
hardtype = round(hardtype);

% Number of half-cycles
ncyc = length(load);

% Convergence tolerance
tol = 1e-7;

% Find the maximum number of points needed
maxlen=0;
for i=1:ncyc
	if length(load{i})>maxlen; maxlen=length(load{i}); end
end
maxlen = maxlen+1;

% Create struct for saving data
hist			= struct;
hist.etot = nan(maxlen,ncyc); % Total strain
hist.p	  = nan(maxlen,ncyc); % Accumulated plastic strain
hist.s    = nan(maxlen,ncyc); % Total stress
hist.a    = nan(maxlen,ncyc); % Total back-stress
hist.R    = nan(maxlen,ncyc); % Yield surface expansion/contraction

% Save a backup for exceptions
hist_backup = hist;

% Initial state
alpha_0		= zeros(1,n_alpha+1);
R_0				= 0;
e_0				= 0;
pstrain_0 = 0;
p_0				= 0;
first			= 1; % Correction factor for initial loading cycle

% Main loop
for idx=1:ncyc
	% Strain points in half-cycle
	eprobe = load{idx}';

	% Positive or negative cycle
	sgn = round(2*(mod(idx,2)-0.5));
	
	% Elastic limit
	eran = (2-first)*(s_y+R_0)/E; % Elastic strain range
	elim = e_0+sgn*eran;					% Elastoplastic transition point

	% Points in elastic/plastic ranges
	eprobe_elast = eprobe(sgn*eprobe<=sgn*elim);
	eprobe_plast = eprobe(sgn*eprobe>sgn*elim);
	nplast = length(eprobe_plast);
	nelast = length(eprobe_elast);
	ntot   = nplast + nelast;

  % If response is purely elastic, no plastic evolution
  if nplast == 0
    pguess  = [];
    delta_p = [];
    eguess  = [];
    alpha   = [];
		R				= [];
	else
		% Map out 100 log spaced points on plastic strain axis, used to 
		% evaluate the resulting total strain in the first iteration. pguess is
		% difference in plastic strain, can be positive and negative.
		pguess = sgn*logspace(-6,log10(abs(eprobe_plast(end)-elim)),100)';
    	
		% Convergence loop
		nits		= 0;
    p       = zeros(length(eprobe),1);
    delta_p = zeros(length(eprobe_plast),1);
    eguess  = zeros(length(eprobe_plast),1);
    alpha   = zeros(length(pguess),n_alpha+1);
    R       = R_0;
		while nits<1000
			nits=nits+1;
	
			% Eq. plastic strain increments 
			delta_p = abs(pguess);
			
			% Calculate backstress
      alpha   = zeros(length(pguess),n_alpha+1);
			for i=1:n_alpha
				alpha(:,i+1)=c(i)/gamma(i)*pguess./delta_p.*(1- ...
				exp(-delta_p*gamma(i)))+exp(-delta_p*gamma(i))*alpha_0(i+1);
			end
			alpha(:,1) = sum(alpha(:,2:end),2); 
			
			% Calculate R and dR := ∂R/∂p = R_dot/p_dot
			if hardtype == 1 % Voce hardening
				R = Q*(1-exp(-delta_p*b))+exp(-delta_p*b)*R_0;
        dR = b*(Q-R);
	
			elseif hardtype == 2 % Power hardening
				R = Q*(p_0+delta_p).^b;
        dR = R*b./((R/Q).^(1/b)); 
	
      else  % Modified power hardening
				v  = 10*(abs(Q*b/sum(c)))^(1/(1-b));
       	R  = Q*(p_0+delta_p + v).^b - Q*v^b;
       	dR = (R+Q*v^b)*b./((R/Q+v^b).^(1/b)); 
			end
	
			% Corresponding elastic and total strain
      eguess = (alpha(:,1)+sgn*(s_y+R))/E;
			guess = (eguess+pguess+pstrain_0);
			  	
			% Evaluate and update plastic strain guess
      if nits == 1
       	% The total strain from the 100 log spaced pstrain points is
				% used to linearly interpolate an approximate position for each
				% data point.
	
				% interp1 needs non-nan, sorted data points.
				nanmask = ~isnan(guess);
				x = guess(nanmask); y = pguess(nanmask);
        [~,donde] = unique(x(:,1),'sorted');
        x = x(donde); y = y(donde);
            	
				% If data points become nan for some reason, throw exception
        if length(x)<2
					hist = hist_backup; 
		     	return;
        end
				
				% Linearly interpolate good pstrain starting points
        pguess = interp1(x,y,eprobe_plast,"linear","extrap");
	
			else
				% The starting points from the previous iteration are used as
				% Taylor expansion points of the equation:
				% e_tot = e_p + (alpha+sgn*(R+s_y))/E
				% A 2nd order expansion is used for backstress, and first order 
				% expansion for R.
			
        % If guess error is below tolerance, accept, otherwise update
        dev = abs(guess-eprobe_plast);
        if all(dev<tol); break; end
            	
				% Calculate first- and second order alpha derivatives
        gradvec = zeros(length(pguess),2);
		    for i=1:n_alpha 
					% Calculate ∂alpha/∂p = alpha_dot/p_dot
         	gradvec(:,1) = gradvec(:,1)+(c(i)*sgn-gamma(i).*alpha(:,i+1));
					% Calculate ∂^2alpha/∂p^2
          gradvec(:,2) = gradvec(:,2)+(-gamma(i))* ...
												 (c(i)*sgn-gamma(i).*alpha(:,i+1));
		    end
            	
				% Collect terms of the 2nd order polynomial: 
				% 0 = C1*p^2 + C2*p + C3
        C1 = 1/2*gradvec(:,2)/E;
        C2 = sgn+(-delta_p.*gradvec(:,2)+gradvec(:,1)+dR*sgn)/E;
        C3 = pstrain_0-eprobe_plast+ (alpha(:,1)-gradvec(:,1).*delta_p +...
                	1/2*delta_p.^2.*gradvec(:,2)+sgn*(s_y+R-dR.*delta_p))/E;
            	
				% Calculate discriminant
        D = C2.^2-4*C1.*C3;
            	
				% If D is negative or imaginary, throw exception
        if ~isempty(D(D<0)) | sum(abs(imag(D)))>0
					hist = hist_backup; 
		     	return;
        end
				
				% Calculate the two solutions to the equation
        solvec = sgn*[(-C2+sqrt(D))./(2*C1), (-C2-sqrt(D))./(2*C1)];
            	
				% Find the solution closest to the initial guess
        difvec = abs(solvec-pguess);
        [~,loc]=min(difvec,[],2);
	
				% Update guess
        pguess = solvec(nplast*(loc-1)+[1:1:nplast]');  
      end
		end
		% If convergence loop exceeded limit, throw exception
		if nits==1000
			hist = hist_backup; 
			return; 
		end
	end

	% Save plastic strain
  pstrain = [zeros(nelast,1); pguess] + pstrain_0;
	p = [zeros(nelast,1); delta_p] + p_0;
	
	% Save data in struct
	hist.etot(1:ntot,idx) = eprobe;
	hist.p(1:ntot,idx)		= p;
	hist.s(1:ntot,idx)		= [eprobe_elast-pstrain_0; eguess]*E;
	if nplast>0
		hist.a(1:ntot,idx)		= [ones(nelast,1)*alpha_0(1); alpha(:,1)];
		hist.R(1:ntot,idx)		= [ones(nelast,1)*R_0; R];
	else
		hist.a(1:ntot,idx)		= [ones(nelast,1)*alpha_0(1)];
		hist.R(1:ntot,idx)		= [ones(nelast,1)*R_0];
	end

	% Update old state
	if ~isempty(pguess) 
		alpha_0		= alpha(end,:);
		R_0				= R(end);
		pstrain_0 = pstrain(end);
		p_0				= p(end); 
	end
	e_0		= eprobe(end);
	first = 0;
end
% Save history to .mat file
if ~isempty(outname)
	save(['OUTDATA/',outname,'.mat'],'-struct','hist');
end
end