function crack
%% CRACK Problem
%
% crack solves the Poisson equation $-\Delta u =f$ in $\Omega$ and $u =
% g_D$ on $\partial \Omega$ in a crack domain
% $\Omega=\{|x|+|y|<1\}\backslash \{0<= x <=1, y=0\}$
%  using adaptive finite element method (AFEM). We choose f=1 and g_D such
%  that the exact solution is $u = r^{\beta}\sin(\beta\theta)-0.25r^2,
%  \beta = 1/2$ in the polar coordinate.
%
% EXAMPLE
%
%    crack 
%
% See also  crack_performance, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all
%% Parameters
maxN = 2e3;     theta = 0.5;    maxIt = 50; 
N = zeros(maxIt,1);     errL2 = zeros(maxIt,1);     errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
elem = label(node,elem);                        % label the mesh
bdFlag = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition
showmesh(node,elem);                            % plot mesh
findelem(node,elem);                            % plot element indices
findnode(node,1:5);                             % plot node indices
findedge(node,[1 5; 5 6],1);                    % plot the crack edge
text(node(6,1)+0.05,node(6,2)+0.075,int2str(6), ...
     'FontSize',14,'FontWeight','bold');
%% 
% node 1 and node 6 are the same point (1,0)

%%  Get a fine mesh by uniform bisection
for k = 1:2
    [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end
clf; showmesh(node,elem);
%% Set up PDE data
pde.f = @f;
pde.g_D = @exactu;
pde.Du=[];
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    u = Poisson(node,elem,pde,bdFlag);
    % Plot mesh and solution
    figure(1);  showresult(node,elem,u,[-7,12]);    
    % Step 2: ESTIMATE
    %eta = estimaterecovery(node,elem,u);            % recovery type
    eta = estimateresidual(node,elem,u,pde);    % residual type
    % Record error and number of vertices
    errL2(k) = getL2error(node,elem,@exactu,u);
    errH1(k) = getH1error(node,elem,@Du,u);
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
end

%% Plot convergence rates
N= N(1:k);  errH1 = errH1(1:k);     errL2 = errL2(1:k); 
figure;
showrate2(N,errH1,3,'-*','||Du-Du_h||',...
          N,errL2,3,'k-+','||u-u_h||');
%%
% Using AFEM, we obtain optimal rate of convergence the error in the energy
% norm ($N^{-0.5}$) and in the $L^2$ norm ($N^{-1}$).
end % End of function CRACK


%%  Data of CRACK
function z = f(p)   % load data (right hand side function)
z = ones(size(p,1),1);
end

function z = exactu(p) % exact solution
r = sqrt(sum(p.^2,2));
z = sqrt(0.5*(r-p(:,1)))-0.25*r.^2;
end

function z = Du(p) % derivative of the exact solution
r = sqrt(sum(p.^2,2));
z(:,1) = (p(:,1)./r-1)./sqrt(8*(r-p(:,1)))-0.5*p(:,1); % u_x
z(:,2) = p(:,2)./r./sqrt(8*(r-p(:,1)))-0.5*p(:,2);     % u_y
end