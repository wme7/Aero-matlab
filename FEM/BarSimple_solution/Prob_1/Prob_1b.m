% 1D Clamped Bar with a Point Load
% Isoparametric Formulation Implementation
% clear memory
clear all;
close all;
clc;
 
% E: modulus of elasticity
% A: area of cross section
% L: length of bar
E=2.E11; A=12.5e-4; L=[0.75 0.75];  
 
% numberElements: number of elements
numberElements=2; 
% numberNodes: number of nodes
numberNodes=3;
% generation of coordinates and connectivities
elementNodes=[1 2;2 3];
nodeCoordinates=[0 0.75 1.5];
% for structure:
   % displacements: displacement vector
   % force : force vector
   % stiffness: stiffness matrix
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes); 
 
% computation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  detJacobian=L(e)/2;
  invJacobian=1/detJacobian;
  ngp = 3;
  [w,xi]=gauss1d(ngp);
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
          stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof)=force(elementDof)+...
          +20000*shape'*detJacobian*w(ip);
  end
end 
% boundary conditions and solution
% prescribed dofs
prescribedDof=[1 3];
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactionsPretty(displacements,stiffness, ...
    numberNodes,prescribedDof,force)
%exact solution
h1=figure();
x_ext=0:0.01:1.5;
u_ext=20000*(1.5-x_ext).*x_ext/A/E/2;
plot(x_ext,u_ext,'r-')
hold on;
plot(nodeCoordinates,displacements,'s--')
legend('Exact solution','FEM');
title('Displacement','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'displacement, $\it{u}$'},'interpreter','latex','FontSize',14);
%evaluate axial stress
h2=figure();
%stress/strain recover
fprintf('Axial stress\n')
fprintf('element\tgauss point\taxial stress\n')
ngp = 2;
elementIpCoor=zeros(numberElements*ngp,1);
outputNode=[];
elementAxiStr=zeros(numberElements*ngp,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  detJacobian=L(e)/2;
  invJacobian=1/detJacobian;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(2)));
  outputNode=[outputNode nodeCoordinates(elementDof(1)) nodeCoordinates(elementDof(2))];
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      elementIpCoor(ip+(e-1)*2,1)=xc+xi(ip)*detJacobian;
      elementAxiStr(ip+(e-1)*2,1)=E*B*displacements(elementDof,1);
      fprintf('%2.0f\t%2.0fth ip\t%10.4e\n', e, ip,elementAxiStr(ip+(e-1)*2,1))
  end
end 
s_ext=20000*(1.5-2*x_ext)/A/2;
plot(x_ext,s_ext,'r-')
hold on;
%plot(elementIpCoor,elementAxiStr,'*')
plot(outputNode,elementAxiStr,'--')
legend('Exact','FEM');
title('Stress','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'Stress, $\it{\sigma}$'},'interpreter','latex','FontSize',14);