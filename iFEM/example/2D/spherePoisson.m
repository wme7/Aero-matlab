function spherePoisson
%% Spheric Poisson Problem
%
% spherePoisson solves the Poisson equation $-\Delta_S u =f$ on $S$, where
% $\Delta _S$ is the surface Laplace-Beltrami operator. We choose $S$ as
% the unit sphere and $f=2x$. The exact solution is $u=x$.
%
% See also crack, Lshape, cube
%
% <a href="matlab:ifemdoc spherePoisson">iFEMdoc spherePoisson</a>
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


maxIt = 5; err = zeros(maxIt,1); N = zeros(maxIt,1);
%% Generate initial mesh
node = [1,0,0; 0,1,0; -1,0,0; 0,-1,0; 0,0,1; 0,0,-1];
elem = [6,1,2; 6,2,3; 6,3,4; 6,4,1; 5,1,4; 5,3,4; 5,3,2; 5,2,1];
showmesh(node,elem); findnode3(node);
for i = 1:2
    [node,elem] = uniformrefine(node,elem);
end

%% 
% project the mesh to the sphere
r = sqrt(node(:,1).^2 + node(:,2).^2 + node(:,3).^2);
node = node./[r r r];
showmesh(node,elem);

%% Finite Element Methods
figure(1); 
set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.7,0.4]);
for k=1:maxIt
    [u,A] = surfacePoisson(node,elem,[],@f,[],[]);
    % Plot mesh and solution
    subplot(1,2,1);
    showmesh(node,elem,[130,28]); pause(0.1);
    subplot(1,2,2); 
    showsolution(node,elem,u,[130,28]); colorbar;
    % Record error and number of vertices    
    uI = exactu(node);
    err(k) = sqrt((u-uI)'*A*(u-uI));
    N(k) = size(node,1);
    % Uniform refinement and Project to the sphere
    [node,elem] = uniformrefine(node,elem);
    r = sqrt(node(:,1).^2+node(:,2).^2+node(:,3).^2);
    node = node./[r r r];
end

%% Plot Error
figure(2)
c = err(ceil(k/2))/N(ceil(k/2))^(-1);
loglog(N,err,'-*',N,c*N.^(-1),'r--','linewidth',2);
title('Superconvergence', 'FontSize', 14);
legend('||Du_I-Du_h||','N^{-1}','LOCATION','Best')
xlabel('Number of vertices N'); 
ylabel('Error');
%%
% The error of the nodal interpolation and the finite element approximation
% in the energy norm is only first order, i.e. $\|u-u_I\|_A \leq Ch,
% \|u-u_h\|_A \leq Ch$. But these two discrete functions are superclose
% $\|u_I - u_h\|_A\leq Ch^2$.

%%
% c = {'DOF','Error'};
% makeHtmlTable([N err],[],[],c,[],[15 3]);
end %End of SPEREPOISSON

%% Data of SPEREPOISSON
function z = f(p) % load data (right hand side function)
z = 2*p(:,1);
end

function z = exactu(p)   % exact solution
z = p(:,1);
end