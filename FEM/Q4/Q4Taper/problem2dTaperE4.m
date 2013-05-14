% Clamped Taper Plate with Vertical Traction
% Q4 Isoparametric Formulation Implementation
% 4 elements
% clear memory
clear all; 
 
% materials
E  = 3*1e7;     poisson = 0.30;  thickness = 1;

% matrix D
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
 
% Brute force preprocessing
% numberElements: number of elements
numberElements=4; 
% numberNodes: number of nodes
numberNodes=9;
% coordinates and connectivities
elementNodes=[1 2 5 4; 2 3 6 5; 4 5 8 7;5 6 9 8];
nodeCoordinates=[0,0; 1,0.25; 2,0.5; 0,0.5; 1,0.625; 2,0.75; 0,1; 1,1; 2,1];
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-o');

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% boundary conditions 
prescribedDof=[1 2 7 8 13 14]';
% force vector 
force=zeros(GDof,1);
force(14)=-10; force(16) = -20; force(18)=-10;

% calculation of the system stiffness matrix
stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,D,thickness);

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements
outputDisplacements(displacements, numberNodes, GDof);

