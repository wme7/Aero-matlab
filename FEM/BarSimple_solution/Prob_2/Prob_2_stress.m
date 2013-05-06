% Isoparametric Formulation Implementation
% clear memory
clear all
close all
clc
% E: modulus of elasticity
% A: area of cross section
% L: length of bar
E=2e5; A=12.5e-4; L=1.5;
s_exact=@(x)1/2*(80000*x)*x/(A);
hold on;
ezplot(s_exact,[0 1.5])
title('Exact Solution v.s. FEM','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('Nodal displacement','interpreter','latex');
for i=1:4   %For NEL = 1, 2, 4, 8 
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
NNOD=2;
numberElements=2^(i-1); 
% numberNodes: number of nodes
numberNodes=numberElements+1;
% generation of coordinates and connectivities           
elementNodes=[1:numberNodes-1;2:numberNodes]';                              
nodeCoordinates=linspace(0,L,numberNodes); 
%Generate element length vector
Le=ones(1,2^(i-1))*L/numberElements; 
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
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  ngp = 2;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(2)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof)=force(elementDof)+...
          -80000*(xc+detJacobian*xi(ip))*shape'*detJacobian*w(ip);
  end
end 
 
% boundary conditions and solution
% prescribed dofs
prescribedDof=[numberNodes];
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactionsPretty(displacements,stiffness, ...
    numberNodes,prescribedDof,force)
%stress/strain recover
fprintf('Axial stress\n')
fprintf('element\tgauss point\taxial stress\n')
ngp = 2;
elementIpCoor=zeros(numberElements*ngp,1);
elementNodeCoor=zeros(numberElements*NNOD,1);
elementAxiStr=zeros(numberElements*ngp,1);
elementNodeStr=zeros(numberElements*NNOD,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(2)));
  Coefficient=zeros(ngp,ngp);
  for ip=1:ngp
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      elementIpCoor(ip+(e-1)*ngp,1)=xc+xi(ip)*detJacobian;
      elementAxiStr(ip+(e-1)*ngp,1)=E*B*displacements(elementDof,1);
      [shape,naturalDerivatives]=shapeFunctionL2(1/xi(ip));
      Coefficient(:,ip)=shape;
      fprintf('%2.0f\t%2.0fth ip\t%10.4e\n', e, ip,...
          elementAxiStr(ip+(e-1)*ngp,1))
  end
  elementNodeCoor(1+(e-1)*ngp:2+(e-1)*ngp)=[nodeCoordinates(elementDof(1))...
      nodeCoordinates(elementDof(2))];
  sigma=elementAxiStr(1+(e-1)*ngp:2+(e-1)*ngp,1);
  elementNodeStr(1+(e-1)*ngp:2+(e-1)*ngp,1)=Coefficient*sigma;  
end 

%post process
switch numberElements
    case 1
        str='rs--';
    case 2
        str='gs--';
    case 4
        str='ks--';
    case 8
        str='ms--';
end
switch numberElements
    case 1
        strg='r*';
    case 2
        strg='g*';
    case 4
        strg='k*';
    case 8
        strg='m*';
end
plot(elementNodeCoor,elementNodeStr,str)
hold on;
%plot(elementIpCoor,elementAxiStr,strg)
end
legend('Exact solution','NEL=1','NEL=2','NEL=4','NEL=8');
title('Stress','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'Stress, $\it{\sigma}$'},'interpreter','latex','FontSize',14);
box on;
