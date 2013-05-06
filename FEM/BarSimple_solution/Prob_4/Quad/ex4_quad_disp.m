% Isoparametric Formulation Implementation
% clear memory
%clear all
close all
clc
% E: modulus of elasticity
% A: area of cross section
% L: length of bar
E=10E4; A=1; c=1; l=1; L=2*l;

%u_exact = @(x) (36*log(x/2)-4*x-12*heaviside(x-5)*log(x/5)+8)/E;
u_exact = @(x) c/(A*E)*(l^2*x-x.^3/6);
du_exact = @(x) c/(A*E)*(l^2-x.^2/2);

hold on;
ezplot(u_exact,[0 2])
title('Exact Solution v.s. FEM','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('Nodal displacement','interpreter','latex');

for i=1:3   %For NEL = 1, 2, 4, 8 
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
numberElements=2^(i-1); 
% numberNodes: number of nodes
numberNodes=2*numberElements+1;
%A=zeros(1,numberElements);
% generation of coordinates and connectivities           
NNOD=3;                            
nodeCoordinates=linspace(0,L,numberNodes);
%Generate element length vector
for k=1:numberElements
    Le(k)=L/numberElements;
    elementNodes(k,:)=[(k-1)*2+1 (k-1)*2+2 (k-1)*2+3];
    %A(k)=(nodeCoordinates(2*k-1)+Le(k)/2)*2;
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
  ngp = 2;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(end)));
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof)=force(elementDof)+...
          c*xc*shape'*detJacobian*w(ip);
  end
  if(nodeCoordinates(elementDof(end))==2)
     x=(2-xc)/detJacobian;
     [s,n]=shapeFunctionL3(x);
     force(elementDof)= force(elementDof)-(c*l^2/A)*s';
  end
  if(nodeCoordinates(elementDof(1))<2&&...
          nodeCoordinates(elementDof(3))>2)
     x=(2-xc)/detJacobian;
     [s,n]=shapeFunctionL3(x);
     force(elementDof)= 0;
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
ngp = 2;
elementIpCoor=zeros(numberElements*ngp,1);
elementNodeCoor=zeros(numberElements*NNOD,1);
elementAxiStr=zeros(numberElements*ngp,1);
elementNodeStr=zeros(numberElements*NNOD,1);
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  elementDof=elementNodes(e,:);
  detJacobian=Le(e)/2;
  invJacobian=1/detJacobian;
  [w,xi]=gauss1d(ngp);
  xc=0.5*(nodeCoordinates(elementDof(1))+nodeCoordinates(elementDof(3)));
  Coefficient=zeros(NNOD,ngp);
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL3(xi(ip));
      B=naturalDerivatives*invJacobian;
      elementIpCoor(ip+(e-1)*ngp,1)=xc+xi(ip)*detJacobian;
      elementAxiStr(ip+(e-1)*ngp,1)=E*B*displacements(elementDof,1);
      fprintf('%2.0f\t%2.0fth ip\t%10.4e\n', e, ip,elementAxiStr(ip+(e-1)*ngp,1))
  end
  Xi=[1/xi(1) 0 1/xi(2)];
  for iNode=1:NNOD
      [shape,naturalDerivatives]=shapeFunctionL2(Xi(iNode));
      Coefficient(iNode,:)=shape;
      elementNodeCoor(iNode+(e-1)*NNOD)=nodeCoordinates(elementDof(iNode));
  end
  sigma=elementAxiStr(1+(e-1)*ngp:ngp+(e-1)*ngp);
  elementNodeStr(1+(e-1)*NNOD:NNOD+(e-1)*NNOD)=Coefficient*sigma;
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

if numberElements > 1
nodes = elementNodeCoor(1:2:end);
stress = elementNodeStr(1:2:end);
else
nodes = elementNodeCoor;
stress = elementNodeStr;    
end

% Save data for analysis
h(i) = Le(1);
L2_quad(i)= sqrt(sum((u_exact(nodeCoordinates') - displacements).^2));
en_quad(i)= sqrt(sum((du_exact(elementNodeCoor) - (elementNodeStr/(A*E))).^2));
%en_quad(i)= sqrt(sum((du_exact(nodes) - (stress/(A*E))).^2));

% Plot actual displacements
plot(nodeCoordinates,displacements,str)
hold on;
end
legend('Exact solution','NEL=1','NEL=2','NEL=4','NEL=8','interpreter','latex');
box on;

%% Plot Errors
figure(2)
loglog(h,L2_quad,'-sb'); hold on; 
loglog(h,en_quad,'--*b'); hold off;
legend('L_2','en','interpreter','latex',4);
xlabel('log_{10}h','interpreter','tex');
box on;