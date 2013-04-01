% Isoparametric Formulation Implementation
% clear memory
clear all

%% Physical Parameters
% numberElements: number of elements
numberElements = 6; 

% Initial Area
A0 = 12.5e-4; % m^2

% Modulus of elasticity
E = 2e11; % N/m^2

% Length of bar
Ltotal = 1.5; % m 

% Build nodes coordinates and middle points,
nodesPerElement = 3;
delta_x = Ltotal/numberElements;
i = 0:numberElements; 
x = delta_x/2*i'; 
x_mid = (x(3)-x(1))/2 + delta_x*i';

% Evaluate area the center of each element.
% A: area of cross section
A = A0;  %*(1+x_mid/Ltotal); % constant cross section

% Length of each element
L = x(2:end)-x(1:end-1);  
 
%% Build Elements
% numberNodes: number of nodes
numberNodes = numberElements + 1;

% generation of coordinates and connectivities
elementNodes = zeros(numberElements,2); 
elementNodes(1,:) = [1,3];    
for i = 2:numberElements        
    elementNodes(i,:) = elementNodes(i-1,:) + 2;    
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
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof) = force(elementDof)+...
          -80000*(xc+detJacobian*xi(ip))*shape'*detJacobian*w(ip);
  end
end  
 
%% BCs and solution
% prescribed dofs
prescribedDof = find(nodeCoordinates == 1.5); % fixed at node x = 1.5 m
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof,force)

%% Stress
stress = zeros(numberElements,1);
for e=1:numberElements;    
    stress(e) = E*B*displacements(elementNodes(e,:));
end

%% Exact solutions
x = 0:0.1:1.5;
displacements_exact=(90000*x-30000*x.^2)/(E*A);
stress_exact=-((60000*x)-90000)/A;

%% plot figures

% Displacements
subplot(1,2,1); hold on;
plot(nodeCoordinates,displacements,'-.sb')
plot(x,displacements_exact,'-r'); hold off;

% Stress
subplot(1,2,2); hold on;
stress = [stress;stress(end)]; 
stairs(nodeCoordinates,stress,'-.sb')
plot(x,stress_exact,'-r'); hold off;