function crack_withcoarsen
% CRACK solves Poisson equation in a crack domain with AFEM.
% For large number of unknowns, comment out the graphic output showmesh.
%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.
%--------------------------------------------------------------------------

close all; clear all;
figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.7,0.85]);
%---------------------- Parameters ----------------------------------------
theta=0.5; maxN = 4e3; maxIt = 50; N = zeros(maxIt,1);
errL2 = zeros(maxIt,1); errH1 = zeros(maxIt,1);
%---------------------- Generate initial mesh------------------------------
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];    % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];        % elements
elem = label(node,elem);          % label the mesh by the longest edge rule
bdEdge = uint8([1 0 1; 1 0 0; 1 0 0; 1 1 0]);    % boundary edges
N0 = size(node,1);                          % record for coarsen
%---------------------- uniform bisection ---------------------------------
for k=1:3
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end
subplot(2,2,1); showmesh(node,elem); pause(0.1)
%---------------------- Adaptive Finite Element Method --------------------          
for k = 1:maxIt
    %------------------- Step 1: Solve ------------------------------------
    u = Poisson(node,elem,bdEdge,@f,@g_D,@g_N);     % use direct solver
    %------------------- Plot mesh and solution ---------------------------
    subplot(2,2,4); showsolution(node,elem,u,[-10,20]);
    subplot(2,2,2); showsolution(node,elem,u,[0,90]); colorbar; pause(0.1);
    %------------------- Step 2: Estimate ---------------------------------
    eta = estimateW21(node,elem,u);            % recovery type
%    eta = estimateresidual(node,elem,u,@f);    % residual type
    %------------------- record error and number of DOF -------------------
    errL2(k) = computeL2error(node,elem,@exactu,u);
    errH1(k) = computeH1error(node,elem,@u_x,@u_y,u);
    %------------------- Step 3: Mark -------------------------------------
    markedElem = mark(elem,eta,theta);
    %------------------- Step 4.1: Refine ---------------------------------
    [node,elem,bdEdge,HB,tree] = bisect(node,elem,markedElem,bdEdge);
    subplot(2,2,1); showmesh(node,elem); pause(0.1)
    %------------------- Step 4.2: Coarsen --------------------------------
    % Not necessary for stationary problem
    eta = eleminterpolate(eta,tree);
    markedElem = mark(elem,eta,0.25*theta,'COARSEN');
    [node,elem,bdEdge] = coarsen(node,elem,markedElem,N0,bdEdge);
    subplot(2,2,1); showmesh(node,elem); pause(0.1)
    N(k) = size(elem,1);
    if (N(k)>maxN), break; end        
end
%---------------------- Plot convergence rates ----------------------------
N= N(1:k); errH1 = errH1(1:k); errL2 = errL2(1:k);
subplot(2,2,3);
loglog(N,errH1,'-*',N,N.^(-1/2),'r--','linewidth',2); 
hold on
loglog(N,errL2,'k-+',N,N.^(-1),'m--','linewidth',2); 
axis tight;
title('H1 error and L2 error ', 'FontSize', 14);
legend('||Du-Du_h||','N^{-0.5}','||u-u_h||','N^{-1}','LOCATION','Best')
xlabel('number of elements NT');
ylabel('error');
end  
%---------------------- End of function CRACK -----------------------------

%---------------------- Sub functions called by CRACK ---------------------
function z = f(p) % load data (right hand side function)
z = ones(size(p,1),1);
end
%--------------------------------------------------------------------------
function z = g_D(p) % Dirichlet boundary condition
z = exactu(p);
end
%--------------------------------------------------------------------------
function z = g_N(p) % Neumann boundary condition
z = zeros(size(p,1),1);
end
%--------------------------------------------------------------------------
function z = exactu(p) % exact solution
r = sqrt(sum(p.^2,2));
z = sqrt(0.5*(r-p(:,1)))-0.25*r.^2;
end
%--------------------------------------------------------------------------
function z = u_x(p) % x-derivative of the exact solution
r = sqrt(sum(p.^2,2));
z = (p(:,1)./r-1)./sqrt(8*(r-p(:,1)))-0.5*p(:,1);
end
%--------------------------------------------------------------------------
function z = u_y(p) % y-derivative of the exact solution
r = sqrt(sum(p.^2,2));
z = p(:,2)./r./sqrt(8*(r-p(:,1)))-0.5*p(:,2);
end
%--------------------------------------------------------------------------