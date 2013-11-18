%% JUMPMG2 Diffusion equation with jump coefficients in three dimensions 
%
% The purpose of this subroutine is to test the performance of V-cycle
% multigrid for solving the linear system of algebraic equations arising
% from linear finite element discretization of the elliptic partial
% differential equation with jump coefficients.
%
% $-\nabla \cdot (\omega\nabla u) = f $  in  $\Omega=(-1,1)^3$ 
% $u = 1$ on $x==1$ and $u=0$ on $x==-1$, 
% $\omega\nabla u \cdot n = 0 on other faces.
%
% The diffusion coefficent $\omega$ is piecewise constant with large jump.
%
% jumpmgdata1
%* $\omega(x) = 1/\epsilon$ if $x\in (0, 1)^3$ and 
%* $\omega = 1$ otherwise.  
%
% jumpmgdata2
%* $\omega(x) = 1$ if $x\in (-0.5, 0)^3$ or $x\in (0,0.5)^3$ and 
%* $\omega = \epsilon$ otherwise.  
%
% See also jumpMG1  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;
%% Input and parameters
maxN = 3e4;     theta = 0.5;      maxIt = 20;      gflag = 1;    
global epsilon
epsilon = 1e-4;

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');

%%  Get a fine mesh by uniform bisection
for k = 1:3
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
end

%% Get the data of the pde
% pde = jumpmgdata1;
pde = jumpmgdata2;

%% Adaptive Finite Element Method
N = zeros(maxIt,1); itStep = zeros(maxIt,1);  cost = zeros(maxIt,1);
for k = 1:maxIt
    option.solver = 'V'; % Vcycle multigrid
    [u,A,eqn,info] = Poisson3(node,elem,pde,HB,bdFlag,option);
    if gflag == 1
    	showresult3(node,elem,u,'~(x<=0 & y<=0)'); 
    end
    itStep(k) = info.itStep;
    cost(k) = info.solverTime;
    N(k) = length(u); 
    if (N(k)>maxN), break; end        
    eta = estimateresidual3(node,elem,u,pde);
    markedElem = mark(elem,eta,theta); 
    [node,elem,HB,bdFlag] = bisect3(node,elem,markedElem,HB,bdFlag);
end

%% Show results
N = N(1:k); itStep = itStep(1:k); cost = cost(1:k);
figure; showrate(N,cost,5)
display(['  DOF     ' '  Iteration   ' '     CPU (second)']); 
display(num2str([N itStep cost]));
%%
% Multigrid alone cannot be robust for jump coefficients problems. The
% iteration as well as cpu time increase more than linearly.