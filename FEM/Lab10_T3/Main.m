%Finite Element Method 101-2
%National Taiwan University
%2D Elasticity problem

%%clear memory
close all; clc; clear all;
format long;
%% Load Mesh for Lab 10
[nodeCoordinates,elementNodes]=Mesh_Lab10('T3');
% node coordinates are given in mm
NodePerElement=3;
numberNodes=length(nodeCoordinates);
numberElements=length(elementNodes);

%% Import BCs

% Essential BC's
GDof = 2*numberNodes;
prescribedDof = [1,2,3,4];

% Natural BC's
force = zeros(GDof,1);
force(end) = -10; % 10 [kN]

%% Import material and section properties
E = 3E7; % [GPa]
poisson = 0.3; %[-]
thickness = 1; % [mm]

%% Evalute force vector
%force=formForceVectorT3(GDof,naturalBCs,surfaceOrientation,...
%    elementNodes,nodeCoordinates,P,thickness);

%% Construct Stiffness matrix for T3 element
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,D,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);
scaleFactor=1.E5;
drawingMesh(nodeCoordinates+scaleFactor*[displacements(1:2:2*numberNodes) ...
    displacements(2:2:2*numberNodes)],elementNodes,'T3','r--');

% Computes elements stresses
for e=1:numberElements                           
  numNodePerElement = length(elementNodes(e,:));
  numEDOF = 2*numNodePerElement;
  elementDof=zeros(1,numEDOF);
  for i = 1:numNodePerElement
      elementDof(2*i-1)=2*elementNodes(e,i)-1;
      elementDof(2*i)=2*elementNodes(e,i);   
  end
  
  %  B matrix
  x1 = nodeCoordinates(elementNodes(e,1),1);
  y1 = nodeCoordinates(elementNodes(e,1),2);
  x2 = nodeCoordinates(elementNodes(e,2),1);
  y2 = nodeCoordinates(elementNodes(e,2),2);
  x3 = nodeCoordinates(elementNodes(e,3),1);
  y3 = nodeCoordinates(elementNodes(e,3),2);
  A = 1/2*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
  B = 1/(2*A).*[y2-y3 0 y3-y1 0 y1-y2 0;
                        0 x3-x2 0 x1-x3 0 x2-x1;
                        x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
    
  stress=D*B*displacements(elementDof);
  vonmises=sqrt(0.5*((stress(1)-(stress(2)))^2+(stress(2))^2+(stress(1))^2+6*(stress(3))^2));
  fprintf('\n Stress in element % u \n',e)
  fprintf('Sigma_xx : %0 .6f \n',stress(1))
  fprintf('Sigma_yy : %0 .6f \n',stress(2))
  fprintf('Sigma_xy : %0 .6f \n',stress(3))
  fprintf('Vonmises : %0 .6f \n',vonmises)
end 