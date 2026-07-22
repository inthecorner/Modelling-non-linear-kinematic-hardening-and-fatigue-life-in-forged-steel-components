clear all
close all

% ==================================
%% DEFINE OPTIMISATION PROBLEM
% ==================================
% Output file name
outname = 'optimisation';

% (INPUT) Name and path
path  = 'TESTDATA_PROCESSED/';
tests = {
				 'example_data'
				 };

% (INPUT) Set number of backstresses
n_alpha = 3;

% Isotropic hardening model
% Options:
%		1 - Voce (traditional)
%		2 - Hollomon (not recommended)
%		3 - Modified Hollomon (recommended)
hardtype = 3;

% (INPUT) Define isotropic bounds
%         s_y		Q				b
lb_iso = [300;	-1000;	0];
ub_iso = [500;	0;		1000];

% (INPUT) Define kinematic bounds
%         c_1		y_1			c_2		y_2		c_3		y_3
lb_kin = [1e3;	50;	    1e5;	1e2;	100;	1];
ub_kin = [2e5;	1e3;		6e6;	2e4;	2e4;	1];

% Collect bounds in total array
lb = [lb_kin(1:n_alpha*2); lb_iso];
ub = [ub_kin(1:n_alpha*2); ub_iso];

% (INPUT) Define starting point
%         s_y			Q		b	
X0_iso = [300;	-150;	0.15];
%       c_1	   y_1	c_2		  y_2		c_3		y_3
X0_kin=[45000; 300; 150000; 1200; 2500; 1];

% Collect
X0 = [X0_kin; X0_iso];

% ==================================
%% CONSTRUCT LOAD CELL
% ==================================
len = length(tests);
data = struct;
for j=1:len
  % Load data
  num = ['data',int2str(j)];
	data.(num) = load([path,tests{j},'/struct.mat']);

  % Choose to optimise on only a subset of cycles (works bad)
  data.(num).subset_idx = 1:1:length(data.(num).proc.cyc);
  data.(num).subset = data.(num).proc.cyc(data.(num).subset_idx);

  % Make a vector with cycle numbers without skipping
  data.(num).cycs = 1:1:data.(num).subset(end);
	
  % Initialise the load cell
  data.(num).load = cell(length(data.(num).cycs),1);
    
  % Pack with measured strain points at selected cycles. At intermediate
  % cycles, skip straight to final point
	for i=2:length(data.(num).cycs)
  	if any(i==data.(num).subset)
	  	idx = find(data.(num).proc.cyc==i);
			nanmask = ~isnan(data.(num).proc.strain(:,idx));
			data.(num).load{i-1}(end+1) = data.(num).proc.point(1,idx);
			data.(num).load{i}	= data.(num).proc.strain(nanmask,idx)';
  	else
			data.(num).load{i-1}(end+1) = data.(num).proc.point_avg(1,...
                                        	abs(mod(i,2)-2));
  	end
	end
end

% Scale parameters to be in the range -100 to 100
paramscale = max(abs([ub;lb]),[],1); paramscale(paramscale==0) = 1;
paramscale = 1./paramscale*100;

% Define objective function as sum of squares 
fun = @(x)R2(x./paramscale,data,n_alpha,hardtype);

% ==================================
%% RUN OPTIMISATION
% ==================================

% Initialsise
X = [];
output = [];

% Run optimisation
tStart = tic;
[X,output]= FminconMain(fun,X0.*paramscale,lb.*paramscale,ub.*paramscale);
tEnd = toc(tStart);

% Scael parameters back to true ranges
X = X./paramscale;

% Load data into struct
out = struct;
out.X0        = X0;
out.X         = X;
out.score     = output;
out.time      = tEnd;

% Save struct
save(['OUTDATA/',outname,'.mat'],'-struct','out');

% =========================
%% FUNCTIONS 
% =========================
function [f]= R2(X,data,n_alpha,hardtype)
  % Number of tests
  len = length(fieldnames(data));

	% Hardening model for Chaboche model
	tmp=reshape(X(1:2*n_alpha),2,[]);
	c=tmp(1,:); gamma=tmp(2,:);
	
	% Isotropic hardening
	s_y = X(2*n_alpha+1); % Yield stress
	Q   = X(2*n_alpha+2); % Saturation limit
	b   = X(2*n_alpha+3); % Saturation rate
	v   = [];             % Shift
	
	% Loop over each experiment
	f = 0;
	for i=1:len
	  num    = ['data',int2str(i)];		% Number of current experiment
		E      = data.(num).proc.Eavg;	% Young's modulus in current data
    load   = data.(num).load;				% Strain points
    
		% Run simulation
		try
			hist=total_form(s_y,E,n_alpha,c,gamma,Q,b,load',hardtype,[]);
  	catch
			f = nan;
			return
		end
		
		% Get experimental stress
		stress = data.(num).proc.stress(:,data.(num).subset_idx);
		
		% Determine rows and columns of stress data to compare
		rows = min(size(hist.s,1),size(stress,1));
		cols = data.(num).proc.cyc(data.(num).subset_idx);

		% Find difference between simulated and experiment stress
		residuals = stress(1:rows,2:end)-hist.s(1:rows,cols(2:end));
		residuals_nonan = residuals(~isnan(residuals));

		% Size of dataset
	  sz = numel(residuals_nonan);
		
		% Weight factor from strain amplitude
		scale  = data.(num).proc.eamp*100;

		% Calculate objective function score and add to running sum
	  if sz > 0 
			f = f + scale*sum(sum(residuals_nonan.^2))/sz; 
    else 
			f = nan;
    end
	end
	% Normalise by amount of tests
	f = f/len;
end

function [x,output,history]=FminconMain(fun,x0,lb,ub)
%% History
% history.it(i):  iteration number  fir ith iteration
% history.fval(i): objective volume for ith iteraton
% history.xe(i,:) design variables for ith iteration
history.it = [];  history.fval = [];  history.xe = [];

%% options for fmincon
%% options for fmincon
options =optimoptions('fmincon', ...
											'OutputFcn',@outfun, ...
											'Algorithm','sqp',...											
											'Display','iter', ...
											'UseParallel',false, ...
											'MaxIterations',300, ...
											'MaxFunctionEvaluations',4000,...											
											'StepTolerance',1e-9, ...
    									'FunctionTolerance',1e-9, ...
											'SpecifyObjectiveGradient',false);
%% Perform optimization using fmincon 
[x,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

%% Function to store optimization history
    function stop = outfun(xe,optimValues,state)
        % output function which is used to store data for each step.
        stop = false;
        switch state
            case 'iter'
                history.fval = [history.fval; optimValues.fval];
                history.it = [history.it; optimValues.iteration];
                history.xe = [history.xe;xe];
            otherwise
        end
    end
end