% A Thin Plate Subjected to Uniform Traction
% T3 Implementation
% 2 elements
% clear memory
clear all; 
clc;
close all;
% materials
E  = 30e6;     poisson = 0.30;  thickness = 1;

% matrix D
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
 
% trivial preprocessing
% numberElements: number of elements
chartLen=zeros(1,4);
answerSet=zeros(5,4);
format long;
for i=1:8
    h=[20/2^(i-1) 10/2^(i-1)];
    xx=0:20/2^(i-1):20;
    yy=0:10/2^(i-1):10;
    [NodeX,NodeY]=meshgrid(xx,yy);
    nx=length(xx);
    ny=length(yy);
    numberElements=2*(length(xx)-1)*(length(yy)-1);
    % numberNodes: number of nodes
    numberNodes=nx*ny;
    % coordinates and connectivities
    elementNodes=[];
    for n=1:nx-1
        for m=1:ny-1
            elementNodes=[elementNodes; (n-1)*nx+m,n*nx+m,n*nx+m+1;...
                (n-1)*nx+m,n*nx+m+1,(n-1)*nx+m+1];
        end
    end
    nodeCoordinates=[NodeX(:) NodeY(:)];
    % GDof: global number of degrees of freedom
    GDof=2*numberNodes; 

    % boundary conditions 
    prescribedDof=(1:2*nx)';
    naturalBCs=(1:2:2*(nx-1))+2*(nx-2)*(ny-1);
    surfaceOrientation=[2];
    P=[1000,0];
    
    % force vector 
    force=formForceVectorT3(GDof,naturalBCs,surfaceOrientation,...
        elementNodes,nodeCoordinates,P,thickness);

    % calculation of the system stiffness matrix
    stiffness=formStiffness2D(GDof,numberElements,...
        elementNodes,numberNodes,nodeCoordinates,D,thickness);

    % solution
    displacements=solution(GDof,prescribedDof,stiffness,force);

    % output displacements
    %outputDisplacements(displacements, numberNodes, GDof);

    %outputStress(displacements,numberElements,...
    %    elementNodes,nodeCoordinates,D)
    results=outputResult(displacements,numberElements,...
        elementNodes,nodeCoordinates,D,h);
    
    
    chartLen(1,i)=h(1);
    answerSet(1:2,i)=displacements(end-1:end);
    answerSet(3:end,i)=results;
end
    %post processing
    figure();
    plot(chartLen,answerSet(1,:),'b*')
    title('x displacement','interpreter','latex','FontSize',18);
    xlabel('element size,$\it{h}$(in)','interpreter','latex','FontSize',14);
    ylabel({'$\it{U_{x}}$,in'},'interpreter','latex','FontSize',14);
    
    figure();
    plot(chartLen,answerSet(2,:),'b*')
    title('y displacement','interpreter','latex','FontSize',18);
    xlabel('element size,$\it{h}$(in)','interpreter','latex','FontSize',14);
    ylabel({'$\it{U_{y}}$,in'},'interpreter','latex','FontSize',14);
    
    figure();
    plot(chartLen,answerSet(3,:),'b*')
    title('Stress','interpreter','latex','FontSize',18);
    xlabel('element size,$\it{h}$(in)','interpreter','latex','FontSize',14);
    ylabel({'$\it{\sigma_{xx}}$,psi'},'interpreter','latex','FontSize',14);
    
    figure();
    plot(chartLen,answerSet(4,:),'b*')
    title('Stress','interpreter','latex','FontSize',18);
    xlabel('element size,$\it{h}$(in)','interpreter','latex','FontSize',14);
    ylabel({'$\it{\sigma_{yy}}$,psi'},'interpreter','latex','FontSize',14);
    
    figure();
    plot(chartLen,answerSet(5,:),'b*')
    title('Stress','interpreter','latex','FontSize',18);
    xlabel('element size,$\it{h}$(in)','interpreter','latex','FontSize',14);
    ylabel({'$\it{\sigma_{xy}}$,psi'},'interpreter','latex','FontSize',14);