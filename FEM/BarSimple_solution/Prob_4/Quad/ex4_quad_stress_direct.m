% Isoparametric Formulation Implementation
% clear memory
clear all
close all
clc
% E: modulus of elasticity
% A: area of cross section
% L: length of bar
E=8; L=4;
u_exact=@(x) (56-8*(x-2)-24*heaviside(x-5))/2/x;
hold on;
ezplot(u_exact,[2 6])
title('Exact Solution v.s. FEM','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('Axial stress, $\it{\sigma}_{x}$','interpreter','latex','FontSize',12);

for i=1:4   %For NEL = 1, 2, 4, 8 
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
numberElements=2^(i-1); 
% numberNodes: number of nodes
numberNodes=2*numberElements+1;
A=zeros(1,numberElements);
% generation of coordinates and connectivities           
NNOD=3;                            
nodeCoordinates=linspace(2,L+2,numberNodes);
%Generate element length vector
for i=1:numberElements
    Le(i)=L/numberElements;
    elementNodes(i,:)=[(i-1)*2+1 (i-1)*2+2 (i-1)*2+3];
    A(i)=(nodeCoordinates(2*i-1)+Le(i)/2)*2;
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
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(end)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A(e);
      force(elementDof)=force(elementDof)+...
          8*shape'*detJacobian*w(ip);
  end
  if(nodeCoordinates(elementDof(end))==5)
     x=(5-xc)/detJacobian;
     [s,n]=shapeFunctionL3(x);
     force(elementDof)= force(elementDof)+...
         24*s';
  end
  if(nodeCoordinates(elementDof(1))<5&&...
          nodeCoordinates(elementDof(3))>5)
     x=(5-xc)/detJacobian;
     [s,n]=shapeFunctionL3(x);
     force(elementDof)= force(elementDof)+...
         24*s';
  end
end 
% boundary conditions and solution
% prescribed dofs
prescribedDof=[1];
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactionsPretty(displacements,stiffness, ...
    numberNodes,prescribedDof,force)
fprintf('Axial stress\n')
fprintf('element\t\taxial stress\n')
%Stress and strain recovery
ngp = 3;
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
        str='r*--';
    case 2
        str='g*--';
    case 4
        str='k*--';
    case 8
        str='m*--';
end

plot(elementNodeCoor,elementNodeStr,str)
hold on;
end
legend('Exact solution','NEL=1','NEL=2','NEL=4','NEL=8','interpreter','latex');
