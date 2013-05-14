% Clamped Taper Plate with Vertical Traction
% Q4 Isoparametric Formulation Implementation
% 1 element
% clear memory
clear all; 
 
% materials
E  = 3e7;     poisson = 0.30;  thickness = 1;

% matrix D
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

% coordinates and connectivities 
x=[0 2 2 0];
y=[0 0.5 1 1];
% Number of partitions 
for power = [1,2,3,4];
% Elements and Connectivities
[numberElements,elementNodes,numberNodes,nodeCoordinates,h,l,partitions] = Q4mesh(x,y,power);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% calculation of the system stiffness matrix
stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,D,thickness);

% boundary conditions 
%prescribedDof=[1 2 3 4]';
prescribedNodes=find(nodeCoordinates(:,1)==0);
xDofs = 2*prescribedNodes-1;
yDofs = 2*prescribedNodes;
prescribedDof = [xDofs;yDofs];

% force vector 
force=zeros(GDof,1);
naturalBCs=find(nodeCoordinates(:,2)==1);
xDofs = 2*naturalBCs-1;
yDofs = 2*naturalBCs;

% traction function
traction = @(x) -20-5*x; %[N/m]
% traction vector
t = traction(nodeCoordinates(naturalBCs,1));

% build force vector
for e = 1:partitions
    force(yDofs(e)) = force(yDofs(e)) + (3*t(e)+t(e+1))*thickness*l/4;
    force(yDofs(e+1)) = force(yDofs(e+1)) + (t(e)+3*t(e+1))*thickness*l/4;
end

% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

% output displacements
outputDisplacements(displacements, numberNodes, GDof);

% Identify data to plot
plotNodes = [partitions+1, numberNodes]'; % nodes (2,1) & (2,0.5)
xDofs = 2*plotNodes-1;
yDofs = 2*plotNodes;

% save data
dx(power) = h;
uy21(power) = displacements(yDofs(1));
uy205(power) = displacements(yDofs(2));

% Draw mesh for each case
figure(1)
subplot(2,2,power); drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
end % end power

% plot
figure(2)
subplot(1,2,1); plot(dx,uy21,'-*'); 
title('U\_y at (2,1)','interpreter','latex');
ylabel('U\_y \textit{(m)}','interpreter','latex');
xlabel('element size \textit{h (m)}','interpreter','latex');
subplot(1,2,2); plot(dx,uy205,'-*'); 
title('U\_y at (2,0.5)','interpreter','latex');
ylabel('U\_y \textit{(m)}','interpreter','latex');
xlabel('element size \textit{h (m)}','interpreter','latex');
