function [sigma,u] = dmg(M,B,C,f,g,elem,option)
%% DMG Distributive Multigrid
%
% [sigma,u] = dmg(M,B,C,f,g,elem) solves saddle point problem discretized
% from mixed FEM.
% 
%      |M   B'| |sigma|  = |f|
%      |B  -C | |u|      = |g|
%
% Use  |D   B'| as the preconditioner in gmres and compute the inverse by
%      |B  -C | 
%
% the factorization
%
% |D B'| |I Dinv*B'| = |D    0       |
% |B -C| |0      -I|   |B B*Dinv*B'+C|
%
% The elliptic matrix B*Dinv*B'+C on the dual grid is inverted by one 
% Vcycle using mg.
%
% It is faster than uzawapcg since no need to compute Minv by an inner pcg
% iteration.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tic;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; tol = option.tol; maxIt = option.solvermaxit; 

%% Set up auxiliary matrices
N = size(M,1); Ng = length(g);
Bt = B';
% DM = spdiags(diag(M),0,N,N);
if isempty(C)
   C = sparse(Ng,Ng);
end
DMinv = spdiags(1./diag(M),0,N,N);
Abar = B*DMinv*Bt + C;

%% Set up matrices in each level
option.solver = 'NO';
[x,info,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Abar,g,elem,option); %#ok<ASGLU>

%% Form big matrix equation
bigA = [M Bt; B -C];
bigF = [f; g];

%% Preconditioned GMRES
option.solver = 'Vcycle';
option.setupflag = false;
option.solvermaxIt = 1;
option.printlevel = 0;
fprintf('Distributive MG #dof: %8.0u   ',Ndof);
x = gmres(bigA,bigF,20,tol,maxIt,@DBC,[],x0);
% x = PFGMRES(bigA,bigF,x0,maxIt,20,tol,@DBC,1);
sigma = x(1:N); u = x(N+1:end);

    function s = DBC(r)
        rf = r(1:N); rg = r(N+1:end);
        ds = DMinv*rf;
%         du = Abar\(rg - B*ds);
%         du = amg(Abar,rg - B*ds,option);
        du = mg(Abar,rg - B*ds,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        s = [ds + DMinv*(Bt*du); -du];
    end
toc;
end