%Finite Element Method 101-2
%National Taiwan University
%2D Elasticity problem

%%clear memory
close all; clc; clear all;

%% Load Mesh for Lab 10
[nodeCoordinates,elementNodes]=Mesh_Lab10('Q8');
% node coordinates are given in mm
NodePerElement=8;
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);

%% Import BCs

% Essential BC's
GDof = 2*numberNodes;
prescribedDof = [1,2,3,4,5,6];

% Natural BC's
force = zeros(GDof,1);
force(end) = -10; % 10 [kN]

%% Import material and section properties
E = 3E7; % [GPa]
poisson = 0.3; %[-]
thickness = 1; % [mm]

%% Evalute force vector
%force=formForceVectorQ4(GDof,naturalBCs,surfaceOrientation,...
%    elementNodes,nodeCoordinates,thickness);

%% Construct Stiffness matrix for Q8 element
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,D,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);
scaleFactor=1.E5;
drawingMesh(nodeCoordinates+scaleFactor*[displacements(1:2:2*numberNodes) ...
    displacements(2:2:2*numberNodes)],elementNodes,'Q4','r--');

%% B matrix & strain
% 3 by 3 quadrature
[gaussWeights,gaussLocations]=gauss2d('3x3');

for e=1:numberElements                           
  numNodePerElement = length(elementNodes(e,:));
  numEDOF = 2*numNodePerElement;
  elementDof=zeros(1,numEDOF);
  for i = 1:numNodePerElement
      elementDof(2*i-1)=2*elementNodes(e,i)-1;
      elementDof(2*i)=2*elementNodes(e,i);   
  end
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)    
      GaussPoint=gaussLocations(q,:);
      xi=GaussPoint(1);
      eta=GaussPoint(2);
    
% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ8(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(elementNodes(e,:),:),naturalDerivatives);
    
%  B matrix
    B=zeros(3,numEDOF);
    B(1,1:2:numEDOF)       = XYderivatives(:,1)';        
    B(2,2:2:numEDOF)  = XYderivatives(:,2)';
    B(3,1:2:numEDOF)       = XYderivatives(:,2)';
    B(3,2:2:numEDOF)  = XYderivatives(:,1)';
   
  elementNodes(e,:);
  dis(1:2:16)=displacements([2*(elementNodes(e,:))-1],1);
  dis(2:2:16)=displacements([2*(elementNodes(e,:))],1);
  strain=B*dis';
  stress(:,q)=D*strain;

  end  
  for q=1:size(gaussWeights,1)   
      GaussPoint=gaussLocations(q,:);
      xi=1/GaussPoint(1);
      eta=1/GaussPoint(2);
      [shapeFunction,naturalDerivatives]=shapeFunctionQ8(xi,eta);
      stressxx(e,q)=stress(1,1:8)*shapeFunction;
      stressyy(e,q)=stress(2,1:8)*shapeFunction;
      stressxy(e,q)=stress(3,1:8)*shapeFunction;
      vonmises(e,q)=sqrt(0.5*((stressxx(e,q)-stressyy(e,q))^2+(stressyy(e,q))^2+(stressxx(e,q))^2+6*(stressxy(e,q))^2));
      fprintf('\nStress in element %u\n',e)
      fprintf('\nStress in node %u\n',q)
      fprintf('Sigma_xx : %0.6f\n',stressxx(e,q))
      fprintf('Sigma_yy : %0.6f\n',stressyy(e,q))
      fprintf('Sigma_xy : %0.6f\n',stressxy(e,q))
      fprintf('Vonmises : %0.6f\n',vonmises(e,q)) 
  end
    
end