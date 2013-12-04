%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with Modal DG
%
%               du/dt + df/dx = S(u), for x \in [a,b]
%                 where f = f(u): linear
%
%              coded by Manuel Diaz, NTU, 2012.12.05
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: Tim Warburton, Numerical Partial Differential Equations, Lecture
% Notes 3-15. http://www.caam.rice.edu/~timwar/MA578S03/MA578.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic scalar upwind scheme implementation with RK integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; %clc;

%% Parameters
cfl = 0.02; % CFL condition
tEnd = 2*pi; % final time
K = 7; % degree of accuaracy
nE = 20; % number of elements

%% PreProcess
% Define our Flux function
a=-1; flux = @(w) a*w; 
dflux = @(w) a*ones(size(w));

% Build 1d mesh
xgrid = mesh1d([0 2*pi],nE,'Radau',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vadermonde; lR = tool.legRightEnd; lL = tool.legLeftEnd;
M = tool.MassMatrix; invM = tool.invMassMatrix; Dr = tool.CoefDiffMatrix;
S = M*Dr;

% IC
u0 = IC(x,1);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% times steps
dt = cfl*dx/abs(a);
t = 0:dt:tEnd;

% Load Initial conditions
u = u0; ut = V\u;

% lift IDs
if a >= 0
    direction = 'upwind';
    idthis = 1:nE;
    idneig = [nE,1:nE-1];
else
    direction = 'downwind';
    idthis = [2:nE,1];
    idneig = 1:nE;
end

for tSteps = t
    % ut old
    uo = ut;
    
    for rkstage = 3:-1:1
        % Volume term
        volterm = -a*Dr*uo;
        
        % Surf terms
        sR = flux(lR*uo);
        sL = flux(lL*uo);
        
        % Lift all
        switch direction
            case 'upwind'
                lift = -lL'*sL(:,idthis)+lL'*sR(:,idneig);
            case 'downwind'
                lift = -lR'*sL(:,idthis)+lR'*sR(:,idneig);
        end
        
        % The residual
        df = volterm + invM*lift;
        
        % Stage
        uo = ut + (dt/rkstage)*df/J;
    end
    % Update ut
    ut = uo;
    
    % Plot
    %if rem(it,10) == 0
    u = V*ut;
    plot(x,u0,'-+',x,u,'-'); axis(plotrange); grid on;
    xlabel('x'); ylabel('u'); title('DG-FEM')
    drawnow;
    %end
end