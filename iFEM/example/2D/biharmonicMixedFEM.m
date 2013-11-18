function fourthorder
%% 
%     laplace^2 u = f;   [0,1]^2
%               u = g_D;  
%       Du\cdot n = g_N.  
%
%more information , you can check fourorderdata.m,biharmonicP1.m,biharmonicP2.m
%biharmonicP3.m.
% Created by Jie Zhou.
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
close all
%% Parameters
%maxN = 2e3;     theta = 0.5;    
maxIt = 6; 
N = zeros(maxIt,1);     errL2 = zeros(maxIt,1);     errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
node = [0,0; 0,1; 1,0;1,1];        % nodes
elem = [1,3,2; 4,2,3];                 % elements
elem = label(node,elem);               % label the mesh
bdEdge = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition

                           
% showmesh(node,elem);                            % plot mesh                
% findelem(node,elem,'all','index','color','g');  % plot element indices
% findnode(node,'all','index','color','r');       % plot node indices


%%  Get a fine mesh by uniform bisection
for k = 1:2
     [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);  
%    [node,elem,bdEdge] = bisect(node,elem,'all',bdEdge);
end

showmesh(node,elem);
findnode(node,'all','index','color','r');
findelem(node,elem,'all','index','color','r');

%% Set up PDE data
pde = fourorderdata;


%%       Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
 for k = 1:maxIt
%                   [w,u] = biharmonicP1(node,elem,pde,bdEdge);
                    [w,u] = biharmonicP2(node,elem,pde,bdEdge);
%                   [w,u] = biharmonicP3(node,elem,pde,bdEdge);
               errL2(k) = getL2error(node,elem,pde.exactu,u);
               errH1(k) = getH1error(node,elem,pde.Du,u);
                   N(k) = size(w,1)+size(u,1);
     [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
 end

%  Plot convergence rates
N= N(1:k);  errH1 = errH1(1:k);     errL2 = errL2(1:k); 
figure;
showrate2(N,errH1,3,'-*','||Du-Du_h||',...
          N,errL2,3,'k-+','||u-u_h||');



