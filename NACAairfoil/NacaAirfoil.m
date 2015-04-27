function [y,x] = NacaAirfoil(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        NACA 4 digit airfoil generator
%                   coded by Manuel Diaz, NTU, 2015.04.27
%
% Refs:
% [1] Geometry for Aerodynamicists, Appendix A. Accessible online at:
%     http://www.dept.aoe.vt.edu/~mason/Mason_f/CAtxtAppA.pdf 
%
% Notation:
%             _ _   _________________
% thickness(t) |   /                 ----------__________
%             _|_ (______________________________________---------
%                 |                                              |
%                 |-------------------chord(c)-------------------|
% 
% Usage: 
%       [xi,yi] = NACAairfoil('MPXX')               % defaul usage, 
%       [xi,yi] = NACAairfoil('MPXX',TE)            % minimal usage, 
%       [xi,yi] = NACAairfoil('MPXX',TE,n,xs)       % uncomplete usage,
%       [xi,yi] = NACAairfoil('MPXX',TE,n,xs,path)  % complete usage,
% where:
%           TE  = 'zero' for zero thickness trailing edge,
%               = 'finite' for thickness trailing edge,
%           n   = number of nodes along x-direction,
%           xs  = 'halfcos' for "half cosine x-spacing",
%               = 'uniform' to give "uniform x-spacing",
%           pth = './' for creating airfoil data file in the current folder
%                 (eg. '/naca2412.dat'),
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controling parameters
switch nargin
    case 0
        AF= '0012';	 % designation
        M = AF(1);	 % The maximum value of the mean line as chord/100.
        P = AF(2);   % Is the chordwise position of the maximum chamber as chord/10.
        XX= AF(3:4); % Is the maximum thickness, t/c, in chord precent (x 100%).
        TE= 'zero';  % Trailing edge (TE)
        n = 30;      % x-nodes to describe airfoil
        xs='halfcos';% x-nodes spacing method
        wf= 0;       % writte solution file
    case 1
        AF= varargin{1}; 
        M = AF(1);
        P = AF(2);
        XX= AF(3:4);
        TE= 'zero';
        n = 30;      
        xs= 'halfcos';
        wf= 0;
	case 2
        AF= varargin{1}; 
        M = AF(1);
        P = AF(2);
        XX= AF(3:4);
        TE= varargin{2};
        n = 30;      
        xs= 'halfcos';
        wf= 0;
	case 4
        AF= varargin{1}; 
        M = AF(1);
        P = AF(2);
        XX= AF(3:4);
        TE= varargin{2};
        n = varargin{3};      
        xs= varargin{4};
        wf= 0;
	case 5
        AF= varargin{1}; 
        M = AF(1);
        P = AF(2);
        XX= AF(3:4);
        TE= varargin{2};
        n = varargin{3};      
        xs= varargin{4};
        wf= 1;
        pth=varargin{5};
    otherwise
        error('wrong number of input arguments, check again dummy.')
end

% set airfoil main characteristics
m=str2double(M)/100;   % mean line coord
p=str2double(P)/10;    % maximun chamber position
t=str2double(XX)/100;  % relative thickness (t/c)

% Fixed parameters
a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;
switch TE 
    case 'finite'   % For finite thick TE
        a4=-0.1015;  
    case 'zero'     % For zero thick TE
        a4=-0.1036;
    otherwise
        error('not a valid trailing edge')
end

%% Build Airfoil cross section
% Set spacing for discrete values in x-direction
switch xs 
    case 'halfcos'  % Half cosine based spacing
        beta=linspace(0,pi,n)'; x=(0.5*(1-cos(beta))); 
        header=['NACA',AF,':[',num2str(2*n),'panels,Half cosine x-spacing]'];  
    case 'uniform'  % uniform based spacing
        x=linspace(0,1,n)';
        header=['NACA',AF,':[',num2str(2*n),'panels,Uniform x-spacing]'];
    otherwise
        error('not a valid spacing method')
end

% find points before and after the maximum chamber position 'p'
xc1=x(x<=p);  % for t <= p
xc2=x(x> p);  % for t > p
xc=[xc1;xc2]; % concatenate vectors

% compute thickness function yt(x) 
yt=(t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

% 
if p==0	% symmetric airfoil
    % thickness embelope with camber or mean line.
    xu=x;    yu= yt;    % upper surface
    xl=x;    yl=-yt;    % lower surface
    yc=zeros(size(xc)); % camber line function
else	% non-symmetric airfoil
    % camber lines
    yc1=(m/p^2)*(2*p*xc1-xc1.^2);               % for t <= p
    yc2=(m/(1-p)^2)*((1-2*p)+2*p*xc2-xc2.^2);   % for t > p
    yc=[yc1;yc2];       % camber line function yc(x)

    dyc1_dx=(m/p^2)*(2*p-2*xc1);                % for t <= p
    dyc2_dx=(m/(1-p)^2)*(2*p-2*xc2);            % for t > p
    dyc_dx=[dyc1_dx;dyc2_dx];   % camber line slope dyc/dx
    theta=atan(dyc_dx);	% camber line slope angle 

    xu= x-yt.*sin(theta);	yu=yc+yt.*cos(theta); % upper surface
    xl= x+yt.*sin(theta);	yl=yc-yt.*cos(theta); % lower surface
end

% Plot airfoid embelope
plot(xu,yu,'bo-'); hold on  % upper surface
plot(xl,yl,'ro-'); hold off % lower surface
axis equal

% Output points
x=[flipud(xu) ; xl(2:end)];
y=[flipud(yu) ; yl(2:end)];


%% Write Points to file
if wf == 1
    F1=header;
    F2=num2str([x,y]);
    F=strvcat(F1,F2);
    fileName=[pth,'NACA',AF,'.dat'];
    dlmwrite(fileName,F,'delimiter','')
end

end