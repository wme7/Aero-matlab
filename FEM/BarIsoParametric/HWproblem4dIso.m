% Isoparametric Formulation Implementation
% clear memory
clear all; 

figure(5)

% Displacements
subplot(1,2,1); hold on;
%plot(xx,displacements_exact,'-r'); 

% Stress
subplot(1,2,2); hold on;
%plot(xx,stress_exact,'-r');

%% Physical Parameters
% Initial Area
A0 = 2; % m^2

% Modulus of elasticity
E = 8; % N/m^2

% Length of bar
Ltotal = 4; % m 

%% Build Elements
% number of elements,
%numberElements = 4; 
test = [1 2 4 8];
for numberElements = test

% number of nodes,
numberNodes = 2*numberElements+1;

% node coordinates,
x = linspace(2,6,numberNodes)';

% Generation of coordinates and connectivities
L = zeros(numberElements,1); elementNodes = zeros(numberElements,3); 
L(1) = Ltotal/numberElements; elementNodes(1,:) = [1,2,3];
for i = 2:numberElements        
    L(i) = L(1); elementNodes(i,:) = elementNodes(i-1,:) + 2;    
end
nodeCoordinates = x;

% element center coordinate
x_mid = 0.5*(x(elementNodes(:,3))-x(elementNodes(:,1)))+x(elementNodes(:,1));

% Evaluate area the center of each element.
% A: area of cross section
A = A0*x_mid;

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
  ngp = 3;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(3)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A(e);
      force(elementDof) = force(elementDof)+...
          8*(xc+detJacobian*xi(ip))*shape'*detJacobian*w(ip);
  end
end  
 
%% BCs and solution
% prescribed dofs
prescribedDof = find(nodeCoordinates == 2); % fixed at node x = 2 m
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness, ...
    numberNodes,prescribedDof,force)

%% Stress
elementCoords=zeros(numberElements*ngp,1);
stress=zeros(numberElements*ngp,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:); 
  detJacobian=L(e)/2;
  invJacobian=1/detJacobian;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(3)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip));
      B=naturalDerivatives*invJacobian;
      elementCoords(ip+(e-1)*3,1)=xc+xi(ip)*detJacobian;
      stress(ip+(e-1)*3,1)=E*B*displacements(elementDof,1);
  end
end 

%% plot figures

switch numberElements
    case 1
        line='r*--';
    case 2
        line='g*--';
    case 4
        line='k*--';
    case 8
        line='m*--';
end

% Displacements
subplot(1,2,1); plot(nodeCoordinates,displacements,line)

% Stress
subplot(1,2,2); axis([2,6,0,26]); plot(elementCoords,stress,line)

end %end for
subplot(1,2,1); legend('Exact','NEL=1','NEL=2','NEL=4','NEL=8',4); hold off
subplot(1,2,2); legend('Exact','NEL=1','NEL=2','NEL=4','NEL=8',1); hold off
