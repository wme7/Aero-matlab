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
title('Exact solution v.s. FEM','interpreter','latex');
xlabel('x,\it{m}','interpreter','latex');
ylabel('Nodal displacement,\it{m}','interpreter','latex');

for i=1:3  %For NEL = 1, 2, 4, 8 
fprintf( '\nNumber of elements:%d\n\n',2^(i-1) );
% numberElements: number of elements
numberElements=2^(i-1); 
% numberNodes: number of nodes
numberNodes=numberElements+1;
% generation of coordinates and connectivities           
elementNodes=[1:numberNodes-1;2:numberNodes]';                              
nodeCoordinates=linspace(0,L,numberNodes);
NNOD=2;
%Generate element length vector
Le=ones(1,numberElements)*L/numberElements;
%A=(nodeCoordinates(1:end-1)+Le./2)*2;

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
  if(nodeCoordinates(elementDof(end))==2)
     x=(2-xc)/detJacobian;
     [s,n]=shapeFunctionL2(x);
     force(elementDof)= force(elementDof)-(c*l^2/A)*s';
  elseif(nodeCoordinates(elementDof(1))<2&&...
          nodeCoordinates(elementDof(end))>2)
     x=(2-xc)/detJacobian;
     [s,n]=shapeFunctionL2(x);
     force(elementDof)= 0;
  end
  for ip=1:ngp;
      [shape,naturalDerivatives]=shapeFunctionL2(xi(ip)); 
      B=naturalDerivatives*invJacobian;
      stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+ B'*B*w(ip)*detJacobian*E*A;
      force(elementDof)=force(elementDof)+...
          c*xc*shape'*detJacobian*w(ip);
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

%stress/strain recovery
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

for e = numberElements
    xe = nodeCoordinates(elementNodes(e,:));
    strain_exact(e,1) = quad(du_exact,xe(1),xe(2))/Le(e);
end
stress = elementNodeStr(1:2:end);

% Save data for analysis
h(i) = Le(1);
L2_linear(i)= sqrt(sum((u_exact(nodeCoordinates') - displacements).^2));
en_linear(i)= sqrt(sum((du_exact(elementNodeCoor) - elementNodeStr/(A*E)).^2));

% Matlab Style
%L2_linear(i) = norm((u_exact(nodeCoordinates') - displacements),2);
%en_linear(i) = norm((du_exact(elementNodeCoor) - elementNodeStr/(A*E)),2);

% Plot actual displacements
plot(nodeCoordinates,displacements,str)
hold on;
end
legend('Exact solution','NEL=1','NEL=2','NEL=4','NEL=8','interpreter','latex',2);
box on;

%% Plot Errors
figure(2)
loglog(h,L2_linear,'-sr'); hold on; 
loglog(h,en_linear,'--*r'); hold off;
legend('L_2','en','interpreter','latex',4);
xlabel('log_{10}h','interpreter','tex');
box on;