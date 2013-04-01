function [displacements,stress,x] =  Driver1d(numberElements)
% 1D Simple Bar code
% Modified by Manuel Diaz 2013.03.16
% close all; clear all; clc;

% Define numberof Elements to use,
% numberElements = 6; 

%% Parameters
% Initial Area
A0 = 12.5/(100)^2; % m^2

% Modulus of elasticity
E = 70e9; % N/m^2

% Length of bar
Ltotal = 0.5; % m 

% Build nodes coordinates and middle points,
delta_x = Ltotal/numberElements;
i = 0:numberElements; x = delta_x*i';
x_mid = (x(2:end)-x(1:end-1))/2 + x(1:end-1);

% Evaluate area the center of each element.
% A: area of cross section
A = A0*(1+x_mid/Ltotal);

% Length of each element
L = x(2:end)-x(1:end-1);

%% Main code

% numberNodes: number of nodes
numberNodes = numberElements+1;

% generation of coordinates and connectivities
elementNodes = zeros(numberElements,2); elementNodes(1,:) = [1,2];
    for i = 2:numberElements
        elementNodes(i,:) = elementNodes(i-1,:) + 1;
    end
nodeCoordinates = x;

% for structure:
   % displacements: displacement vector
   % force : force vector
   % stiffness: stiffness matrix
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes); 
% applied load at node 2
force(1)=-5000.0;
 
% computation of the system stiffness matrix
k = zeros(1,numberElements);
    for e=1:numberElements; 
      % elementDof: element degrees of freedom (Dof)
      elementDof=elementNodes(e,:) ;
      k(e)=E*A(e)./L(e) ;
      stiffness(elementDof,elementDof)=...
           stiffness(elementDof,elementDof)+k(e)*[1 -1;-1 1];
    end
% Boundary conditions and solution

% prescribed Dofs
%prescribedDof=[1;4];
prescribedDof = find(nodeCoordinates == 0.5); % fixed at node x = 0.5 m

% solution
GDof=numberNodes;
displacements = solution(GDof,prescribedDof,stiffness,force);

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof,force);

% output stress
stress = zeros(numberElements,1);
for e=1:numberElements;
    stress(e) = E*1./L(e)*[-1 1]*displacements(elementNodes(e,:));
end

return
