function cubePoissonP2
%% CUBEPOISSONP2 solves Poisson equation in a cube using quadratic element.
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Parameters 
maxN = 3e5; maxIt = 3; N = zeros(maxIt,1); 
error = zeros(maxIt,1); cost = zeros(maxIt,1);

%% Generate initial mesh
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
bdFlag = setboundary3(node,elem,'Dirichlet');
N0 = size(node,1);
HB(1:N0,1:3) = repmat((1:N0)',1,3); HB(1:N0,4) = 0;

%%  Get a fine mesh by uniform bisection
for k = 1:2
    [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
end

%% Solve and Refine
for k=1:maxIt
    t = cputime;
    [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
%     [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    %[node,elem] = uniformbisect3(node,elem);    
    %[u,A,b,edge] = Poisson3P2(node,elem,[],@f,@g_D,[]); 
    [u,A,b,edge] = Poisson3P2(node,elem,bdFlag,HB,@f,@g_D,[]); 
    cost(k) = cputime - t;
    N(k) = size(node,1) + size(edge,1);
    uI = exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    error(k) = (u-uI)'*A*(u-uI);
    %error(k) = norm(u(1:N(k))-uI(1:N(k)));
    if (N(k)>maxN), break; end        
end

%% Graphic output
cost = cost(1:k); N = N(1:k); error = error(1:k);
figure(1); 
r1 = showrate(N,cost,2,'-*'); display(r1)
figure(2);
r2 = showrate(N,error,2,'k-+');
legend('||DuI-Du_h||',['N^{' num2str(r2) '}'], 'LOCATION','Best');
end
%---------------------- End of function CUBE ------------------------------

%---------------------- Sub functions called by CUBE ----------------------
function z=f(p) % load data (right hand side function)
x = p(:,1); y = p(:,2); z = p(:,3);
r2 = x.^2 + y.^2 + z.^2;
z = 20*exp(-10*(r2)).*(3-20*r2);
end
%--------------------------------------------------------------------------
function z = g_D(p) % Dirichlet boundary condition
z = exactu(p);
end
%--------------------------------------------------------------------------
function z = exactu(p) % exact solution
normp = sum(p.^2,2);
z = exp(-10*normp);
end