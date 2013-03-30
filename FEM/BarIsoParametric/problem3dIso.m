% Isoparametric Formulation Implementation
% clear memory
clear all

%% Parameters
% numberElements: number of elements
numberElements = 6; 

% Initial Area
A0 = 12.5e-4; % m^2

% Modulus of elasticity
E = 2e11; % N/m^2

% Length of bar
Ltotal = 1.5; % m 

% Build nodes coordinates and middle points,
delta_x = Ltotal/numberElements;
i = 0:numberElements; 
x = delta_x*i';x_mid = (x(2:end)-x(1:end-1))/2 + x(1:end-1);

% Evaluate area the center of each element.
% A: area of cross section
A = A0;  %*(1+x_mid/Ltotal); % constant cross section

% Length of each element
L = x(2:end)-x(1:end-1);  
 
%% Main code 
% numberNodes: number of nodes
numberNodes = numberElements + 1;

% generation of coordinates and connectivities
elementNodes = zeros(numberElements,2); 
elementNodes(1,:) = [1,2];    
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
% computation of the system stiffness matrix and force vector
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  detJacobian=L(e)/2;
  invJacobian=1/detJacobian;
  ngp = 2;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(2)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof) = force(elementDof)+...
          -80000*(xc+detJacobian*xi(ip))*shape'*detJacobian*w(ip);
  end
end  
 
% boundary conditions and solution
% prescribed dofs
prescribedDof = find(nodeCoordinates == 1.5); % fixed at node x = 1.5 m
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof,force)

% output stress
stress = zeros(numberElements,1);
for e=1:numberElements;    
    stress(e) = E*1./L(e)*[-1 1]*displacements(elementNodes(e,:));
end

%% plot figures
% Displacements
subplot(1,2,1); 
plot(nodeCoordinates,displacements,'-.sb')

% Stress
subplot(1,2,2)
stress = [stress;stress(end)];
stairs(nodeCoordinates,stress,'-.sb')