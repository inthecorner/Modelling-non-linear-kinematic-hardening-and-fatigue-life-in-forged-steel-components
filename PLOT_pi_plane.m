clear all
close all

% ===================================
% SETTINGS
% ===================================
% Animate (1) or trace (2)
plotstyle = 1;

% If animating, how fast?
skip = 50;

% If tracing, how many frames?
ntrace = 20;

% If tracing, draw intermediate stress vectors?
interm = 0;

% Draw gridlines?
gridlines = 1;

% Create GIF?
gif = 0;

% Use principal stresses (1) or diagonal values (0)?
prin = 1;

% How many datapoints do we trim off the end?
remove = 0;

% ===================================
% READ DATA
% ===================================
data = load("OUTDATA\incremental_strain.mat");
len = length(data.smax)-remove;
s_y = data.smax(1); % Inital yield stress
outname = data.name;

% ===================================
% MAKE PI-PLANE
% ===================================
[fig,leg] = piplane(s_y,gridlines,'\sigma',1);

% ===================================
% MAKE THE GUTS
% ===================================
% Yield surface circle
theta = linspace(0,2*pi,100);
x = cos(theta); y = sin(theta);

% Movement directions in plot (unit vectors)
dir1 = sqrt(2/3)*[cos(90*pi/180);sin(90*pi/180)];
dir2 = sqrt(2/3)*[cos(330*pi/180);sin(330*pi/180)];
dir3 = sqrt(2/3)*[cos(210*pi/180);sin(210*pi/180)];

% Make initial empty plots
col4 = [0.4660    0.6740    0.1880];
p4 = scatter(0,0,Marker='o',MarkerFaceColor=col4,MarkerEdgeColor='k',...
						 LineWidth=0.6,SizeData=50);
p1 = plot(0,0,'k');
p2 = plot(0,0,LineStyle="-", color='red',Marker='o',MarkerFaceColor='k',...
					MarkerEdgeColor='k',LineWidth=2);
p3 = quiver(0,0,0,0,s_y*0.3,LineWidth=1,color='#80B3FF');
p3.MaxHeadSize = 2;
title('Yield surface movement in the \pi-plane')

% Save colormap for tracing
cmap = colormap("jet");

% Pick frames depending on plot style
if plotstyle == 1
	interval = skip;
else
	interval = floor(len/ntrace);
end

% Animation loop
nits = 0;
for i=1:interval:len
	nits = nits+1;
	
	% Radius of yield surface
	r = sqrt(2/3)*data.smax(i);
	
	% Get data
	if prin==0
		alpha = [data.a11(i), data.a22(i), data.a33(i)];
		sigma = [data.s11(i), data.s22(i), data.s33(i)];
		de    = [data.edot11(i), data.edot22(i), data.edot33(i)];
	elseif prin == 1
		alpha = [data.ap1(i), data.ap2(i), data.ap3(i)];
		sigma = [data.sp1(i), data.sp2(i), data.sp3(i)];
		de    = [data.edotp1(i), data.edotp2(i), data.edotp3(i)];
	end

	% Make pi-plane vectors
	center      = dir1*alpha(1) + dir2*alpha(2) + dir3*alpha(3);
	stresspoint = dir1*sigma(1) + dir2*sigma(2) + dir3*sigma(3);
	straindir   = dir1*de(1)    + dir2*de(2)    + dir3*de(3);
	
	% If tracing, make a plot per frame, otherwise just update existing
	if plotstyle == 2
		p1 = plot(0,0,color=cmap(1+(nits-1)*floor(size(cmap,1)/ntrace),:));
		p2 = plot(0,0,LineStyle="-", color='red',Marker='o', ...
		MarkerFaceColor='k',MarkerEdgeColor='k',LineWidth=1);
		p3 = quiver(0,0,0,0,s_y*0.3,LineWidth=1,color='#80B3FF');
		p3.MaxHeadSize = 2;
		cbar = colorbar;
		p4 = scatter(0,0,Marker='o',MarkerFaceColor=col4,MarkerEdgeColor=col4);
	end
	
	% Make plots
	p1.XData = x*r + center(1); 
	p1.YData = y*r + center(2);
	
	if plotstyle==2 & interm == 0 & 0<i & i<len
  else
		p2.XData = [0,stresspoint(1)];
	  p2.YData = [0,stresspoint(2)];  
  end
  p4.XData = [0,center(1)];
	p4.YData = [0,center(2)];
  p3.XData = stresspoint(1);
	p3.YData = stresspoint(2);
	p3.UData = straindir(1)/norm(straindir);
	p3.VData = straindir(2)/norm(straindir);
	
	drawnow
	pause(0.1)
	if gif==1
		frame = getframe(fig);
		im{nits} = frame2im(frame);
	end
end
legend(leg{1:end},'Yield surface centre','','','Plastic strain direction','interpreter','none',Location='south')
cbar.Label.String='Pseudo-time';cbar.Label.FontSize=12;

% ===================================
%% CREATE GIF
% ===================================
if gif ==1 
	filename = ['OUTDATA/',outname,'.gif']; % Specify the output file name
	nImages = length(im);
	for idx = 1:nImages
		[A,map] = rgb2ind(im{idx},256);
    if idx == 1
			imwrite(A,map,filename,"gif",LoopCount=Inf, DelayTime=0.02);
    else
      imwrite(A,map,filename,"gif",WriteMode="append", DelayTime=0.02);
		end
	end
end

% ===================================
%% HELPER FUNCTIONS
% ===================================
function [fig,leg] = piplane(s_y,gridlines,symb,figscale)
leg = {};	
% ===================================
	% PLOT ARROWS
	% ===================================
	fig = figure(1);
	screenSize = get(0, 'ScreenSize');
	screenW = screenSize(3);
	screenH = screenSize(4);
	set(gcf, 'Position',  [floor(screenW/2)-300, floor(screenH/2)-300, 600, 600])
	set(gcf,'Color','white')
	
	% Parameters
	r = 1.5*s_y;       % axis length
	theta = deg2rad([90 -30 -150]);
	axis off; hold on;
	axis equal; pbaspect([1 1 1]); 
	
	headLen = 0.1 * r;   % length of arrowhead (along shaft)
	headWid = 0.04 * r;   % half-width of arrowhead
	
	% Draw grid lines
	tickFrac = [0 0.25 0.5 0.75]; % fractions of r
	u = [cos(theta); sin(theta)]; % 2x3 matrix, each column is a unit vector
	% Draw grid lines first (light, dashed). For each axis k, draw lines
	% parallel to axis k through ticks on the next axis (to avoid duplicates).
	if gridlines == 1
		gridColor = [0.9 0.9 0.9];
		gridStyle = '-';
		L = 2.4 * r; % extend length of grid lines
		for k = 1:3
    		j = mod(k,3) + 1; % next axis index
    		for t = tickFrac
        		p = t * r * u(:,j);                 % point on axis j
        		Pstart = p - L * u(:,k);
        		Pend   = p + L * u(:,k);
        		plot([Pstart(1), Pend(1)], [Pstart(2), Pend(2)], 'Color', gridColor, 'LineStyle', gridStyle);
				Pstart = -p - L * u(:,k);
        		Pend   = -p + L * u(:,k);
        		plot([Pstart(1), Pend(1)], [Pstart(2), Pend(2)], 'Color', gridColor, 'LineStyle', gridStyle);
				leg = {leg{1:end},'',''};
    		end
		end
	end
	
	% Draw axes as lines + filled triangular heads
	for k = 1:3
    	x = r*cos(theta(k));
    	y = r*sin(theta(k));
    	% Shaft (without quiver head)
    	plot([0, x*(1-headLen/r)], [0, y*(1-headLen/r)], '-k', 'LineWidth', 1.5);

    	% Unit direction and perpendicular
    	u = [cos(theta(k)); sin(theta(k))];
    	perp = [-u(2); u(1)];
	
    	% Triangle points for filled head (tip at end)
    	P1 = [x; y];
    	P2 = P1 - headLen * u + headWid * perp;
    	P3 = P1 - headLen * u - headWid * perp;
	
    	fill([P1(1), P2(1), P3(1)], [P1(2), P2(2), P3(2)], 'k', 'EdgeColor', 'none');
		leg = {leg{1:end},'',''};
	end
	
	% Add ticks and labels
	tickFrac = [0.25 0.5 0.75];
	tickLen = 0.02 * r;
	% Axis labels
	for k = 1:3
    	th = theta(k);
    	for t = tickFrac
        	px = t*r*cos(th);
        	py = t*r*sin(th);
        	dx = -tickLen*sin(th);
        	dy =  tickLen*cos(th);
        	plot([px-dx, px+dx], [py-dy, py+dy], 'k', 'LineWidth', 1);
			leg = {leg{1:end},''};
    	end
    	text(1.08*r*cos(th), 1.08*r*sin(th), [symb,'_',int2str(k)], ...
         	'HorizontalAlignment','center', 'FontSize',12);
	end
	xlim([-r r]*1.2*figscale); ylim([-r*0.7 r]*1.2*figscale);
end