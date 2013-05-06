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
title('Exact solution v.s. FEM','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('Axial stress, $\it{\sigma}_{x}$','interpreter','latex');
for i=1:4  %For NEL = 1, 2, 4, 8 
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
numberElements=2^(i-1); 
% numberNodes: number of nodes
NNOD = 2;
numberNodes=numberElements+1;
% generation of coordinates and connectivities           
elementNodes=[1:numberNodes-1;2:numberNodes]';                              
nodeCoordinates=linspace(2,L+2,numberNodes);
%Generate element length vector
Le=ones(1,numberElements)*L/numberElements;
A=(nodeCoordinates(1:end-1)+Le./2)*2;
%Evaluate area for each cross-section at the center of each element by A=2x; 
% for structure:
   % displacements: displacement vector
   % force : force vector
   % stiffness: stiffness matrix
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes); 
% computation of the system stiffness matrix and force vector
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:);
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  ngp = 2;
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(end)));
  [w,xi]=gauss1d(ngp);
  if(nodeCoordinates(elementDof(1))==5)
     x=(5-xc)/detJacobian;
     [s,n]=shapeFunctionL2(x);
     force(elementDof)= force(elementDof)+...
         24*s';
  elseif(nodeCoordinates(elementDof(1))<5&&...
          nodeCoordinates(elementDof(end))>5)
     x=(5-xc)/detJacobian;
     [s,n]=shapeFunctionL2(x);
     force(elementDof)= force(elementDof)+...
         24*s';
  end
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A(e);
      force(elementDof)=force(elementDof)+...
          8*shape'*detJacobian*w(ip);
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
%stress/strain recover
fprintf('Axial stress\n')
fprintf('element\tgauss point\taxial stress\n')
ngp = 2;
elementNodeCoor=zeros(numberElements*NNOD,1);
elementNodeStr=zeros(numberElements*NNOD,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:) ;
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  % Nodal coordiantes in natural coordinate
  xi=[-1 1];
  for ip=1:NNOD
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      elementNodeStr(ip+(e-1)*NNOD,1)=E*B*displacements(elementDof,1);  
      fprintf('%2.0f\t%2.0fth ip\t%10.4e\n', e, ip,...
          elementNodeStr(ip+(e-1)*NNOD,1))
  end
  elementNodeCoor(1+(e-1)*ngp:2+(e-1)*ngp)=[nodeCoordinates(elementDof(1))...
      nodeCoordinates(elementDof(2))]; 
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
