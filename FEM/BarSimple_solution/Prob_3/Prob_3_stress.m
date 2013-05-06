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
title('Exact solution v.s. FEM','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('Axial stress, Pa','interpreter','latex');

for i=1:4   %For NEL = 1, 2, 4, 8 
NNOD=3;
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
numberElements=2^(i-1); 
% numberNodes: number of nodes
numberNodes=2*numberElements+1;
% generation of coordinates and connectivities           
                            
nodeCoordinates=linspace(0,L,numberNodes); 
%Generate element length vector
for i=1:numberElements
    Le(i)=L/numberElements;
    elementNodes(i,:)=[(i-1)*2+1 (i-1)*2+2 (i-1)*2+3];
end
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
  ngp = 3;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(3)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip)); 
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
fprintf('Axial stress\n')
fprintf('element\t\taxial stress\n')
%Stress and strain recovery
ngp = 2;

elementNodeCoor=zeros(numberElements*NNOD,1);
elementNodeStr=zeros(numberElements*NNOD,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:);
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  xi=[-1 0 1];
  for ip=1:NNOD;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip));
      B=naturalDerivatives*invJacobian;
      elementNodeCoor(ip+(e-1)*NNOD,1)=nodeCoordinates(elementDof(ip));
      elementNodeStr(ip+(e-1)*NNOD,1)=E*B*displacements(elementDof,1);
      fprintf('%2.0f\t%2.0fth node\t%10.4e\n', e, ip,elementNodeStr(ip+(e-1)*NNOD,1))
  end
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