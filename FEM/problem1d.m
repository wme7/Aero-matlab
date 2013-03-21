%function results =  problem1d(numberElements)
% 1D Simple Bar 
% clear memory
clear all

% numberElements: number of elements
numberElements = 6; 

%% Parameters
A0 = 12.5/(100)^2; % m^2
% E: modulus of elasticity
E = 70e9; % N/m^2

% L: length of bar
Ltotal = 0.5; % m 

% node coordinate and middle points.
delta_x = Ltotal/numberElements;
i = 0:numberElements; x = delta_x*i;
x_mid = (x(2:end)-x(1:end-1))/2 + x(1:end-1);

% Evaluate area the center of each element.
% A: area of cross section
A = A0*(1+x_mid/Ltotal);

% length of each element
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
force(2)=5000.0;
 
% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  k(e)=E*A(e)./L(e) ;
  stiffness(elementDof,elementDof)=...
       stiffness(elementDof,elementDof)+k(e)*[1 -1;-1 1];
end
% boundary conditions and solution
% prescribed dofs
prescribedDof=[1;4];
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof,force)
