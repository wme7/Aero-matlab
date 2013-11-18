function Lshape_uniform
%% LSHAPE_UNIFORM Problem
%
% LSHAPE solves the Poisson equation $-\Delta u =f$ in $\Omega$ and $u =
% g_D$ on $\partial \Omega$ in a crack domain $\Omega=(-1,1)^2\backslash
% (0,1)\times (-1,0)$
%  using adaptive finite element method (AFEM). We choose f and g_D such
%  that the exact solution is $u = r^{\beta}\sin(\beta\theta), \beta = 2/3$
%  in the polar coordinate.
%
% To see the improvement using AFEM, run Lshape_uniform to get the
% convergent rate using uniform refinement.
%
% EXAMPLE
%    Lshape 
%
% See also  crack, Lshape_uniform
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
%% Parameters
maxN = 3e4;     maxIt = 5; 
N = zeros(maxIt,1);   energy = zeros(maxIt,1);  uIuhErrH1 = zeros(maxIt,1);

%%  Generate an initial mesh
node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
bdEdge = setboundary(node,elem,'Dirichlet');
elem = fixorientation(node,elem);   % counter-clockwise oritentation
elem = label(node,elem);            % label the mesh by the longest edge rule
showmesh(node,elem);                % plot mesh
findelem(node,elem);                % plot element indices
findnode(node);                     % plot node indices

%%  Get a fine mesh by uniform bisection
for k = 1:3
    [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
end

clf; showmesh(node,elem);
%% Set up PDE data
pde.f = 0;
pde.g_D = @exactu;
pde.Du=[];% used for recoverFlux;
pde.d=[];
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    [u,A] = Poisson(node,elem,pde);
    % Plot mesh and solution
    N(k) = size(node,1);
    if N(k) < 2e3 % show mesh and solution for small size
        figure(1);  showresult(node,elem,u,[-50,12]);    
    end
    % figure(1);  showresult(node,elem,u(1:size(node,1)),[-50,12]);    
    % Record error and number of vertices
    energy(k) = u'*A*u;
    uI = exactu(node);
    %uI = exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    uIuhErrH1(k) = sqrt((uI-u)'*A*(uI-u));
    if (N(k)>maxN), break; end        
    % Step 4: REFINE
    [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
end

%% Plot convergence rates
N = N(1:k); 
uIuhErrH1 = uIuhErrH1(1:k);
energyError = sqrt(energy(1:k)-energy(k));
figure;
showrate2(N,uIuhErrH1,1,'-*','||Du_I-Du_h||',...
          N(1:k-1),energyError(1:k-1),1,'k-+','sqrt{E(u_k)-E(u_i)}');
%%
% In this example, since f=0, the Dirichlet energy of u is $\|u\|_A^2$. By
% the minimization of the Galerkin projection, we compute $\|u-u_i\|_A^2
% \approx \|u_k - u_i\|_A^2 = \|u_k\|_A^2 -\|u_i\|_A^2$.
%
% We also compute the energy norm between the nodal interpolation and the
% finite element approximation. It is shown that $\|u_I-u_h\|_A$ admits
% convergent rate more than optimal one $N^{-1/2}$. This is known as
% superconvergence. For a finite element function v, the squared energy
% norm can be computed as $\|v\|_A^2 = v'*A*v$.
end % End of function LSHAPE


function z = exactu(p) % exact solution
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
z = r.^(2/3).*sin(2*theta/3);
end

