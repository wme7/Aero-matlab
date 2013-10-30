%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Solving 1-D scalar wave equation with CPR/FR
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Simulation Parameters
CFL = 0.9; % CFL condition
tEnd = 0.2; % final time
kDeg = 3; % degree of accuaracy

parameter.cfl = CFL; 
parameter.tEnd = tEnd; 
parameter.kDeg = kDeg; 
parameter.nP = kDeg+1; 

nE = 5; % number of Elements
range = [0,1]; % Domain range

grid.nE = nE;
grid.nP = nP;
grip.NF = nF;

%% LGL quadrature data
%[xi,w] = GaussLegendre(parameter.nP);
[xi,w,V] = GaussLobatto(parameter.nP);

%% Load right-Radau Polynomials
gR = RadauRightP(parameter.nP); % 1-order higher

%% compute gR'(xi) & gL'(xi)
dgR = diff(gR,sym('x'));
dgR = double(subs(dgR,xi));
dgL = flipud(dgR);

%% Build Elements Nodes & Grid
         
%numberElements = 5;
   numberNodes = numberElements*numberSPs; 
  nodedbCoords = linspace(range(1),range(2),numberElements+1); 
   elemen.tSize = (range(2)-range(1))/numberElements;
 element.Center = (nodedbCoords(1:numberElements) + ...
                  nodedbCoords(2:numberElements+1)) / 2;
        metric = element.Size/2;
[scaledSPs,xc] = meshgrid(elementCenter,metric*normSPs);
    nodeCoords = scaledSPs+xc; % x-Grid
  element.Nodes = reshape(1:number.Nodes,numberSPs,numberElements);
element.dbNodes = [element.Nodes(1,:); element.Nodes(numberSPs,:) ];

%% Load IC
