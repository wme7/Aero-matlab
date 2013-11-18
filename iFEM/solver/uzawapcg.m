function [sigma,u] = uzawapcg(M,B,C,f,g,elem,option)
%% UZAWAPCG solves saddle point problem discretized from mixed FEM.
% 
% [u,sigma] = uzawapcg(M,B,f,g,C) solves the saddle point problem
%      |M  B'| |sigma|  = |f|
%      |B -C | |u|      = |g|
%
%
% Example
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tic;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; tol = option.tol; maxIt = option.solvermaxit; 
% additional parameters for interiori pcg 
tolM = 1e-8; pcgmaxIt = 200; 

%% Set up auxiliary matrices
N = size(M,1); Ng = length(g);
Bt = B';
DM = spdiags(diag(M),0,N,N);
if isempty(C)
   C = sparse(Ng,Ng);
end
Abar = B*spdiags(1./diag(M),0,N,N)*Bt + C;

%% Set up matrices in each level
option.solver = 'NO';
[x,info,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Abar,g,elem,option); %#ok<ASGLU>

%% PCG iteration for solving the Schur complement equation
%
%   (B*M^{-1}*B' + C)u = B*M^{-1}f - g
%
% Preconditioner: B*DM^{-1}*B' + C which is an elliptic matrix on the dual
% grid and can be solved by amg or one Vcycle.

b = B*(1./diag(M).*f) - g; % approximated rhs
nb = norm(b);
% initial residual
u = x0(N+1:end);
[tempr,flag] = pcg(M,f-Bt*u,tolM,pcgmaxIt,DM); % M^{-1}(B'*u)
r = B*tempr - g - C*u;

err = zeros(maxIt,2);
err(1,:) = norm(r)/nb; 
k = 1;
option.solver = 'Vcycle';
option.setupflag = false;
option.solvermaxit = 1;
option.printlevel = 0;
while (max(err(k,:)) > tol) && (k <= maxIt)
    % given r, compute Pr    
%     Pr = Abar\r;  
%     Pr = amg(Abar,r,option);
%     Pr = mg(Abar,r,elem);
    Pr = mg(Abar,r,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
    % update tau, beta, and p
    rho = Pr'*r;  % e'*ABA*e approximates e'*A*e
    if k==1
        p = Pr;
    else
        beta = rho/rho_old;
        p = Pr + beta*p;
    end
    % update alpha, u, and r
    [tempp,flag] = pcg(M,Bt*p,tolM,pcgmaxIt,DM); %#ok<*NASGU> % M^{-1}(B'*p)
    Ap = B*tempp + C*p;  % A*p = B*M^{-1}*B'*p + C*p; 
    alpha = rho/(Ap'*p);
    r = r - alpha*Ap;
    u = u + alpha*p;
    rho_old = rho;
    % compute err for the stopping criterion
    k = k + 1;
    err(k,1) = sqrt(abs(rho/(u'*b))); % approximate relative error in energy norm
    err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
end
fprintf('#dof: %8.0u, Uzawa PCG iter: %2.0u, err = %12.8g\n',...
         Ndof, k, max(err(k,:)));

%% Solve sigma     
[sigma,flag] = pcg(M,f-Bt*u,tol,pcgmaxIt,DM);
toc;