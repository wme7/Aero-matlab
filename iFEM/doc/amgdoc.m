%% ALGEBRAIC MULTIGRID METHOD
%
% We consider solving an SPD matrix equation |Ax = b|, where |A| could be
% obtained by the finite element discretization on a unstructured grids. A
% coarsening of the graph of |A| is needed and restriction and prolongation
% can be constructued based on the coarsening. We focus on the classic
% algebraic multigrid based on the strong connectness. A breif introduction
% can be found at <https://e-reports-ext.llnl.gov/pdf/333205.pdf Falgout, RD. An Introduction to Algebraic Multigrid>

%% Coarsening
% See <./coarsenAMGdoc.html coarsenAMGdoc> 

%% Prolongation
%
% The submatrix A_{cf} is used to construct the interpolation of values on
% fine nodes from that of coarse nodes. The weight is normalized to
% preserve the constant.

%% Test: Different meshes
%
% We consider linear finite element discretization of the Poisson equation
% with homongenous Dirichlet boundary condition on different meshes. We
% summarize the results in <./amgdoctest1.html AMG Test I> which shows our
% AMG solver is of complexity O(Nlog(logN)).

%% Test: Different boundary conditions
%
% We consider the linear finite element discretization of Poisson equation
% on the unstructured mesh with Dirichlet or Neumann boundary conditions.
% We summarize the results in <./amgdoctest2.html AMG Test II> which shows
% our AMG solver is robust to different boundary conditions.

%% Test: Different time stepsize
%
% We consider the linear finite element discretization of heat equation on
% the unstructured mesh with Neumann boundary conditions. We test the
% implicit time discretization with various time stepsizes dt = 1/h^4,
% 1/h^2, 1/h, 1. We summarize the results in <./amgdoctest3.html AMG Test III> which shows
% our AMG solver is robust to different time steps. 

%% Test: Three dimensional problems
%
% We consider the linear finite element discretization of Poisson equation
% in three dimensions with Dirichlet boundary conditions and compare
% geometric multigrid and algebraic multigrid in <./amgdoctest4.html AMG
% Test IV>.
%
% The algebraic multigrid is robust to the size of the matrix. The
% iteration steps increases slightly. The preformance is better on
% structure grids than unstructure grids. The sparsity of the coarse grid
% is increased
