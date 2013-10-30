%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Solving 1-D Burgers' equation with CPR/FR
%
%              coded by Manuel Diaz, KU, 2013.08.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

%% Simulation Parameters
parameter.clf = 0.1/100;	% CFL condition
parameter.tEnd = 0.50;      % End time - Parameter part of ICs
% Polynomial degree
parameter.kDeg = 3;         % Polinomial Degree
parameter.numberSPs = kDeg+1;  % Number of solution points in [-1,1]
% Using RK integration time step
parameter.RK_stages = 3;	% Number of RK stages
    
%% Build Lobatto & Radau Polynomials
% use Lobatto points for the solution points in the standar element [-1,1]
% and a left Radau Polynomial in this implementation, 
       normSPs = rootsLobattoPolynomial(numberSPs);
             g = radaurightPolynomial(k_deg+1); % degree k+1.
            dg = diff(g); % derivative, degree k.

%% Build Elements Nodes & Grid
         range = [0,1];
numberElements = 5;
   numberNodes = numberElements*numberSPs; 
  nodedbCoords = linspace(range(1),range(2),numberElements+1); 
   elementSize = (range(2)-range(1))/numberElements;
 elementCenter = (nodedbCoords(1:numberElements) + ...
                  nodedbCoords(2:numberElements+1)) / 2;
        metric = elementSize/2;
[scaledSPs,xc] = meshgrid(elementCenter,metric*normSPs);
    nodeCoords = scaledSPs+xc; % x-Grid
  elementNodes = reshape(1:numberNodes,numberSPs,numberElements);
elementdbNodes = [elementNodes(1,:); elementNodes(numberSPs,:) ];

%% Time Steps
dt = CFL*elementSize;
steps = tEnd/dt;

%% Define IC
% Assume that the initial condition is u0 = 0.5+sin(2*pi*x)
u0 = 0.5 + sin(2*pi*nodeCoords); 

%% Main Program
% load IC
u = u0; plot(nodeCoords,u,'o-');

%% Math Objects
% Load solutions points for standard element [-1,1]
x0 = normSPs;

% Build Lagrange base functions
x = sym('x');
for i=1:k_deg+1
    l(i)=x/x;
    for j=1:k_deg+1
        if(i ~= j) 
            l(i)=l(i)*(x-x0(j))/(x0(i)-x0(j));
        end
    end
end

% Compute derivatives of lagrange base functions
for i=1:k_deg+1
    D(i)=diff(l(i),x);
    for j=1:k_deg+1
        coefd(j,i) = subs(D(i),x0(j));
    end
end
coefd = double(coefd);

% Compute dg(x0) to the left and right
dgl = subs(dg,x0); dgl = double(dgl); dgr = flipud(dgl);

%% Solver Stage
% Using RK3 time integration,
scale = 1/metric;
% 3-stage Runge Kutta
for n=1:steps
    uo = u;
       
    % 1st stage
    res = burgers1d_CPRrhs(u,k_deg,dgl,dgr,coefd,numberElements,elementdbNodes);
    u = uo-dt*res*scale;

    % 2nd Stage
    res = burgers1d_CPRrhs(u,k_deg,dgl,dgr,coefd,numberElements,elementdbNodes); 
    u = 0.75*uo+0.25*(u-dt*res*scale);

    % 3rd stage
    res = burgers1d_CPRrhs(u,k_deg,dgl,dgr,coefd,numberElements,elementdbNodes); 
    u = (uo+2*(u-dt*res*scale))/3;
    
    %Plot the solution and the exact solution
    plot(nodeCoords,u,'o-');
    
    %pause(0.1)
    if rem(n,2) == 0
        drawnow;
    end
end
