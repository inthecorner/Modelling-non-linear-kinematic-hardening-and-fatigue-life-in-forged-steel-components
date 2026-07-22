clear all
close all

% ==================================
%% (INPUT) DEFINE SIMULATION
% ==================================
% Path and test name
path = 'TESTDATA_PROCESSED/';
test = 	'example_data';

% Load data
data = load([path,test,'/struct.mat']);

% Optimised parameters
X=[
114669.439308535
447.905731352547
1960896.69542957
6576.28197660003
1103.37790000000
1
300
-244.586149109640
0.117931731183189
];

% Isotropic hardening model
% Options:
%		1 - Voce (traditional)
%		2 - Hollomon (not recommended)
%		3 - Modified Hollomon (recommended)
hardtype = 3;

% ==================================
%% DATA PROCESSING
% ==================================
% Average Youngs modulus in test
E = mean(data.proc.young(floor(length(data.proc.young)*1/4): ...
                         floor(length(data.proc.young)*4/5)));

% Number of backstresses
n_alpha = (length(X)-3)/2;

% Extract c and gamma parameters as individual arrays
tmp=reshape(X(1:2*n_alpha),2,[]);
c=tmp(1,:); gamma=tmp(2,:);

% Isotropic hardening parameters
s_y = X(2*n_alpha+1); % Yield stress
Q   = X(2*n_alpha+2); % Saturation limit
b   = X(2*n_alpha+3); % Saturation rate

% Find strain amplitude and mean strain
e_amp  = diff(data.proc.point_avg(1,:))/2;
e_mean = sum(data.proc.point_avg(1,:))/2;
if e_mean<0.2*data.proc.eamp
    e_mean = 0;
end

% Create data for abaqus
cols = rgb2hex(orderedcolors("gem"));
p_sample = [0;logspace(-6,log10(10),499)'];
v = 10*(abs(Q*b/sum(c)))^(1/(1-b));
HardData_ABQ = [s_y + Q*(p_sample + v).^b - Q*v^b,p_sample];

% ==================================
%% RUN SIMULATION
% ==================================
% Initial loading
load{1} = linspace(0,e_mean+e_amp*1,100);
% Main loading
for i=2:2:data.proc.cyc(end)
	load{i}   = linspace(e_mean+e_amp,e_mean-e_amp,100);
	load{i+1} = linspace(e_mean-e_amp,e_mean+e_amp,100);
end

% Run simulation
hist = total_form(s_y,E,n_alpha,c,gamma,Q,b,load,hardtype,[]);

% =============================
%% FIGURE 1 - FULL CYCLES 
% =============================
figure(1)
fig=gcf; fig.Position(3:4)=[920,380];
tiledlayout(1,2);
nexttile
hold on; grid on
x1	 = reshape(hist.etot(:,data.proc.cyc)*100,[],1);
y1	 = real(reshape(hist.s(:,data.proc.cyc),[],1));
z1	 = reshape(zeros(size(y1)),[],1);
col1 = reshape(repmat(data.proc.cyc,size(hist.s,1),1),[],1);
surf([x1,x1]', [y1,y1]', [z1,z1]', [col1,col1]',...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 1);
cbar=colorbar;colormap(jet);
cbar.Label.String='Reversal nr.';cbar.Label.FontSize=12;
xlim([e_mean-data.proc.eamp*1.2,e_mean+data.proc.eamp*1.2]*100)
ylim([-800,800])
xlabel([char(949),'_{11} [%]'])
ylabel('\sigma_{11} [MPa]')
title('Uniaxial cyclic curve (simulation)')
nexttile
hold on; grid on
x2   = reshape(data.proc.strain*100,[],1);
y2   = reshape(data.proc.stress,[],1);
z2   = reshape(zeros(size(y2)),[],1);
col2 = reshape(repmat(data.proc.cyc,size(data.proc.strain,1),1),[],1);
surf([x2,x2]', [y2,y2]', [z2,z2]', [col2,col2]',...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 1);
cbar=colorbar;colormap(jet);
cbar.Label.String='Reversal nr.';cbar.Label.FontSize=12;
xlim([e_mean-data.proc.eamp*1.2,e_mean+data.proc.eamp*1.2]*100)
ylim([-800,800])
xlabel([char(949),'_{11} [%]'])
ylabel('\sigma_{11} [MPa]')
title('Uniaxial cyclic curve (experiment)')

% =============================
%% FIGURE 2 - BACKSTRESSES 
% =============================
delta_p = linspace(0,data.proc.eamp,1000);
alpha = zeros(n_alpha+1,length(delta_p));
for i=1:n_alpha
	alpha(i+1,:)=c(i)/gamma(i)*(1-exp(-delta_p*gamma(i)));
end
alpha(1,:) = sum(alpha(2:end,:),1); 

fig=figure(2);
hold on; grid on
labs = cell(n_alpha+1,1);
for i=1:n_alpha
	plot(delta_p*100,alpha(i+1,:),LineWidth=1)
	labs{i}=['\alpha_',num2str(i)];
end
plot(delta_p*100,alpha(1,:),Color='k',LineWidth=1)
labs{end}= '\alpha_{tot}';
legend(labs,location='northwest')
xlabel('Plastic strain [%]')
ylabel('Stress [MPa]')
title('Backstress development')
fig.Position(3:4)=[410,350];
xlim([0,data.proc.eamp*100]);