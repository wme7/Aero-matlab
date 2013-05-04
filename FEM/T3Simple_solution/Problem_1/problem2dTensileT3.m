% A Thin Plate Subjected to Uniform Traction
% T3 Implementation
% 2 elements
% clear memory
clear all; 
clc;
close all;
% materials
E  = 30e6;     poisson = 0.30;  thickness = 1;

% matrix D
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
 
% trivial preprocessing
% numberElements: number of elements
numberElements=2; 
% numberNodes: number of nodes
numberNodes=4;
% coordinates and connectivities
elementNodes=[1 3 2; 1 4 3];
nodeCoordinates=[0, 0; 0, 10; 20, 10; 20, 0];
drawingMesh(nodeCoordinates,elementNodes,'T3','k-o');

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% boundary conditions 
prescribedDof=[1 2 3 4]';
% force vector 
force=zeros(GDof,1);
force(5)=5000; force(7) =5000;

% calculation of the system stiffness matrix
stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,D,thickness);

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements
outputDisplacements(displacements, numberNodes, GDof);

outputStress(displacements,numberElements,...
    elementNodes,nodeCoordinates,D)