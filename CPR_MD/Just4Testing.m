%% Just fot Testing
% Solving a one-dimensional conservation Law.
% coded by Manuel Diaz, KU, 2013.08.29
clear all; clc; %close all;

%% Simulation Parameters
           CFL = 5/100;   % CFL condition
          tEnd = 0.50;     % End time - Parameter part of ICs
        method = 6;        % {6} CPR only
     plot_figs = 1;        % 0: no, 1: yes please!
     write_ans = 0;        % 0: no, 1: yes please!
% Polynomial degree
         k_deg = 3;        % Polinomial Degree
     numberSPs = k_deg+1;  % Solution Points
% Using RK integration time step
     RK_stages = 3;        % Number of RK stages

%% Define our Flux function
     f = @(w) w.^2/2; %a*w; 
% and the Derivate of the flux function
    df = @(w) w; % a*ones(size(w));
    
%% Build Lobatto & Radau Polynomials
% use Lobatto points and a left Radau Polynomial in this implementation,
       normSPs = rootsLobattoPolynomial(numberSPs);
        g_left = radauPolynomial(k_deg);

%% Build Elements Nodes & Grid
         range = [0,1];
numberElements = 10;
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
              
%% Define IC
% assume the initial condition is u0 = 0.5+sin(2*pi*x)
u0 = 0.5 + sin(2*pi*nodeCoords); %plot(nodeCoords,u0)

%% Main Program
% load IC
u = u0;

% compute flux in every element node,
q = f(u);

% Mask nodes for f_{i+1/2}^{+}
pNodes = elementdbNodes(1,2:end);

% Mask nodes for f_{i+1/2}^{-}
nNodes = elementdbNodes(2,1:end-1);

% compute riemann fluxes (using LF flux spliting)
alpha = max( df(u(pNodes)) - df(u(nNodes)) );
Nfluxes = 0.5*(alpha*(u(pNodes)-u(nNodes)) + (q(pNodes) + q(nNodes)));


% Solution points in the standard element [-1,1]
for i=1:k_deg+1
    x0(i)= -cos((i-1)*pi/k_deg);
end

x = sym('x');

% Lagrange polynomials
for i=1:k_deg+1
    l(i)=x/x;
    for j=1:k_deg+1
        if(i ~= j) 
            l(i)=l(i)*(x-x0(j))/(x0(i)-x0(j));
        end
    end
end
