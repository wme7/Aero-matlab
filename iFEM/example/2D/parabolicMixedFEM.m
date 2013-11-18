function [sigma,u,N,err] = parabolicMixedFEM(dtstr,solver)
%% PARABOLICMIXEDFEM solves the parabolic equation using mixed FEM
%
% parabolicMixedFEM(dtstr,solverstr) discretizes the parabolic equation
% using mixed FEM with step size dt and solves the saddle point equation by
% the solver given by solverstr:
%  'direct': direct solver
%  'uzawapcg': PCG solving the Schur complement equation
%  'dmg': multigrid based on a distributive form
%
% Example
%    parabolicMixedFEM('h^2','uzawapcg');
%    parabolicMixedFEM('1','uzawapcg');
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
for i = 1:5
    [node,elem] = uniformrefine(node,elem);
end
bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
pde = mixBCdata;
maxIt = 3; err = zeros(maxIt,2); N = zeros(maxIt,1);
option.solver = 'none';
for i = 1:maxIt
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    h = 1/sqrt(size(elem,1)); %#ok<NASGU>
    dt = eval(dtstr);
    [tempvar1,tempvar2,eqn] = PoissonRT0(node,elem,pde,bdFlag,option);
    v = simplexvolume(node,elem);
    C = spdiags(v/dt,0,size(elem,1),size(elem,1));
    tic;
    switch lower(solver)
        case 'direct'
            x = A\F;
            sigma = x(1:NE);
            u = x(NE+1:end);
        case 'uzawapcg'
        [sigma,u] = uzawapcg(eqn.M,eqn.B,C,eqn.f,eqn.g,elem);
        case 'dmg'
        [sigma,u] = dmg(eqn.M,eqn.B,C,eqn.f,eqn.g,elem);
    end
    toc;
    N(i) = size(u,1);
end