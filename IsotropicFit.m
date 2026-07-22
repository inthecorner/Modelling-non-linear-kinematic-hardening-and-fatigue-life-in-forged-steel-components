clear all
close all

% =======================================
%% DEFINE OPTIMISATION
% =======================================
% Name and path
path = 'TESTDATA_PROCESSED/';
test = 'example_data';

% Load data
data = load([path,test,'/struct.mat']);

% Isotropic hardening model
% Options:
%		1 - Voce (traditional)
%		2 - Hollomon (not recommended)
%		3 - Modified Hollomon (recommended)
hardtype = 3;

% Initial values for Gauss-Newton optimisation
X	 = [370;-80;0.3];	

% c-values of back-stress components [MPa] (for modified Hollomon)
c = [50000,100000];

% =======================================
%% ESTIMATE ACCUMULATED PLASTIC STRAIN
% =======================================
% Make a vector with cycle numbers without skipping
cycs = 1:1:data.proc.cyc(end);

% Assuming small R changes each cycle, make a vector containing plastic
% strain each cycle
pcum     = zeros(1,data.proc.cyc(end));
pointvec = zeros(2,data.proc.cyc(end));

% Insert stress and strain end points at known cycles
pointvec(:,data.proc.cyc) = data.proc.point; 

% Make a mask specifying the empty spaces from cycle skipping
[~,loc]=intersect(cycs,data.proc.cyc);
empt = zeros(1,length(cycs));empt(loc)=1;empt=logical(empt); empt = ~empt;

% Make a mask to differentiate even and odd cycles
odd = zeros(1,length(cycs)); odd(1:2:end) = 1; odd = logical(odd);

% Make array to contain elastic moduli each reversal, fill in where known
Evec = zeros(1,length(cycs));
Evec(~empt)=data.proc.young;

% In missing cycles, use average strain amplitude as end points
pointvec(1,odd&empt)  = data.proc.point_avg(1,1); 
pointvec(1,~odd&empt) = data.proc.point_avg(1,2);

% Stress decreases each cycle, so average cannot be used. Use linear interp
pointvec(2,odd&empt)  = interp1(data.proc.cyc(1,data.proc.sgn==1), ...
                               data.proc.point(2,data.proc.sgn==1), ...
                               cycs(:,odd&empt),'linear');
pointvec(2,~odd&empt) = interp1(data.proc.cyc(1,data.proc.sgn==-1), ...
                               data.proc.point(2,data.proc.sgn==-1), ...
                               cycs(:,~odd&empt),'linear');
% Also use linear interpolation for elastic moduli in missing spots
Evec(empt) = interp1(data.proc.cyc, data.proc.young,cycs(:,empt),'linear');

% Calculate the plastic strain end the end of each reversal 
pvec = pointvec(1,:)-pointvec(2,:)./Evec;

% Make accumulated plastic strain by summing the difference 
pcum = cumsum([0,abs(diff(pvec))]);

% =======================================
%% RUN OPTIMISATION
% =======================================
% Mask to remove potential points with missing s_y
mask2 = ~isnan(data.proc.sy);

% Make x and y points to fit model to
x = pcum(data.proc.cyc(mask2))';
y = data.proc.sy(mask2)';
% Initial loading should not be considered
x = x(2:end); y = y(2:end); 

% Run optimisation
[X,R2] = benjamin(x,y,X,hardtype,c);

disp(['s_y: ',num2str(round(X(1),3)),' MPa'])
disp(['Q  : ',num2str(round(X(2),3)),' MPa'])
disp(['b  : ',num2str(round(X(3),5))])

% Plot to evaluate fit
fig=figure(7);
scatter(x,y,25,'k.')
hold on; grid on
if hardtype==1
    tit='Voce law';
    plot(x,X(1)+X(2)*(1-exp(-X(3)*x)),LineWidth=1)
elseif hardtype==2
    tit='power law';
    plot(x,X(1)+X(2)*x.^X(3),LineWidth=1);
elseif hardtype==3 
		v = 10*(-X(2)*X(3)/sum(c))^(1/(1-X(3)));
    plot(x,X(1)+X(2)*(x+v).^X(3)-X(2)*v^X(3),LineWidth=1)
    tit = 'modified power law';
end
xlim([0,25])
ylim([250,500])
xlabel('Accumulated plastic strain [-]');ylabel('Yield stress [MPa]')
title(['Fit of experimental yield stress using ',tit])
legend('Experimental data','Least squares model')
fig.Position(3:4)=[420,380];
exportgraphics(fig,"OUTDATA/modhollofit.png",Resolution=600)

% =======================================
%% GAUSS-NEWTON SCRIPT
% =======================================
function [X,R2] = benjamin(x,y,X0,hardtype,c)
	% x,y = data points, column vectors	
	% X0  = initial parameter values [sigma_y; Q; b], column vector
	% Settings and loop variables
	tol = 1e-10;	nmax = 1000;  nits = 0;  err = 1;  R2	= 0;
	% Set initial values
	X	= X0;
	
	% Run convergence loop
	while nits<nmax && err > tol
		nits  = nits+1;
		R2_old = R2;
		% Jacobian & residuals
		if hardtype == 1 % Voce
			J = -1*[ones(length(x),1) 1-exp(-X(3)*x) X(2)*x.*exp(-X(3)*x)];
			r = y-(X(1)+X(2)*(1-exp(-X(3)*x)));
		elseif hardtype == 2 % Hollomon
			J = -1*[ones(length(x),1) x.^X(3) X(2)*x.^X(3).*log(x)];
			r = y-(X(1)+X(2)*x.^X(3));
		elseif hardtype == 3 % Modified Hollomon
      v = 10*(-X(2)*X(3)/sum(c))^(1/(1-X(3)));   	
			J = -1*[ones(length(x),1), (x+v).^X(3)-v^X(3), ...
							X(2)*(x+v).^X(3).*log(x+v)-X(2)*v^X(3)*log(v)];
			r = y-(X(1)+X(2)*(x+v).^X(3)-X(2)*v^X(3));
		end
		
		% Search direction
		h = (J'*J)\(-J'*r);
		
		% Soft line search
		alpha   = logspace(-3,1,9);
		Xmat    = X + alpha.*h;
		
		if hardtype == 1		 % Voce
			Rmat = y-(Xmat(1,:)+Xmat(2,:).*(1-exp(-Xmat(3,:).*x)));
		elseif hardtype == 2 % Hollomon
			Rmat = y-(Xmat(1,:)+Xmat(2,:).*x.^Xmat(3,:));
		elseif hardtype == 3 % Modified Hollomon
      Vmat = 10*(abs(Xmat(2,:).*Xmat(3,:)/sum(c))).^(1./(1-Xmat(3,:)));   	
			Rmat = y-(Xmat(1,:)+Xmat(2,:).*(x+Vmat).^Xmat(3,:)- ...
								Xmat(2,:).*Vmat.^Xmat(3,:));
		end
		Rmat = sum(Rmat.^2,1);
		[~,loc] = min(Rmat);
		damp    = alpha(loc);
		
		% Update guess
		X = X + damp*h;
		
		% Evaluate convergence
		R2  = r'*r/length(r);
		err = abs(R2-R2_old)/R2;
	end
	if nits == nmax; disp('Warning: iteration cap reached'); end
	
	% Calculate R^2 coefficient
	SS_res = sum(r.^2,"all",'omitnan');
	SS_tot = sum((y-mean(y,'omitnan')).^2,1,'omitnan');
	R2 = 1-(SS_res./SS_tot);
end