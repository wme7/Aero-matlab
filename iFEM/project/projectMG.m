%% Project: Multigrid Methods
%
% In this project we will learn three ways of implementating multigrid
% methods: from matrix-free version to matrix-only version depending on how
% much information on the grid and PDE is provided.

%% *Multigrid on Uniform Grids for Poisson Equations*
%
% We consider linear finite element or equivalently 5-point stencil
% discretization of the Poisson equation on a uniform grid of [0,1]^2 with
% size h. For simplicity, we assume h = 1/2^L and zero Dirichlet bounary
% condition.
set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.4,0.3]);
[node,elem] = squaremesh([0,1,0,1],0.25);
subplot(1,2,1); showmesh(node,elem);
[node,elem] = uniformrefine(node,elem);
subplot(1,2,2); showmesh(node,elem);

%% Step 1 Gauss-Seidel iteration
%
% * Use two-dimension array to store functions, e.g. u(i,j) and b(i,j). Then
%   realize A*u and G-S iteration by a double for loop over interior nodes.
%   For example, one Gauss-Seidel iteration is 

n = 2^5; h = 1/n; 
u = zeros(n,n); b = ones(n,n)*h^2;
for i = 2:n-1
   for j = 2:n-1
       u(i,j) = (b(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
   end
end

%%
% * Check the convergence of Gauss-Seidel iteration by solving $-\Delta u = 1$
%   in (0,1)^2 with zero Dirichlet boundary condition. Plot the error in
%  suitable norm vs iteration for h = 1/64.
%
% * Change h from 1/4 to 1/128 and compare the iterations to drive the
%  error in a suitable norm below the discretiztion error h^2.
%
% * Choose random initial guess. Plot the error function on the grid for
% the first 3 steps. Check |surf| and |meshgrid| functions in Matlab for
% the plot of a 2-D function.
%
% * (Optional) Code the red-black Gauss-Seidel by choosing the
% <../redblack.pdf red-black>
% ordering of nodes. (Write two loops by checking if |mod(i+j,2)==0|).
%
%% Step 2 Two-level method
% 
% * Figure out the index map between fine grid with size h and coarse grid
% with size 2h. For example, (i,j) in coarse grid is (2*i-1,2*j-1) in the
% fine grid.
%
% * Code the linear prolongation and restriction using the index map. Be
% carefuly on the value on the boundary points.
%
% * Code the two-grid SSC method. On fine grid, apply m times G-S iteration
% and then restrict the residual to the coarse grid. On coarse grid, use
% G-S iteration in Step 1 to solve the equation below the discretization
% error. Then prolongate the correction to the fine grid and apply
% additional m G-S iterations.
%
% * Change h from 1/4 to 1/128 and compare the iterations of two-level
% methods.
%
%% Step 3 V-cycle Multi-grid method
%
%  Choose one of the following approach to implement the MG.
%
% * (Recrusive way) Apply the two-level method to the coarse grid problem
% in Step 2.
%
% * (Non-recrusive way) Follow the description of SSC in lecture notes to
% implement V-cycle.
%
% _For a non-recrusive implementation of multigrid, the cell structure can
% be used to store the residual and correction with different size on
% different levels._
%
% * Change h from 1/4 to 1/128 and check the iterations and cpu time of MG.

%% *Multigrid on Hierarchical Grids*
%
% We consider linear finite element discretization of the Poisson equation
% on grids obtained by uniform refinement of a coarse grid.
[node,elem] = circlemesh(0,0,1,0.25);
subplot(1,2,1); showmesh(node,elem);
[node,elem] = uniformrefine(node,elem);
subplot(1,2,2); showmesh(node,elem);

%% Step 1 Hierarchical Meshes
%
% Generate the initial grid by 
[node,elem] = circlemesh(0,0,1,0.25);
%%
% * Refine the initial mesh J times to get the finest mesh. To get a mesh of
% the disk, the boundary nodes should be projected onto the unit circle.
%
% * Construct |HB| in two ways. Either from the output of |uniformrefine|
% during the refinement or call |uniformcoarsenred| from the finest mesh.
%
%% Step 2 Transfer Operator and Smoothers
%
% * Assemble the stiffness matrix A in the finest mesh.
%
% * Construct prolongation and restriction matrices using |HB|.
%
% * Compute stiffness matrix in each level by triple product.
%
% * Store the smoother tril(A) and triu(A) in each level.
%
%% Step 3 V-cycle Multigrid
%
% * Follow the lecture notes <http://math.uci.edu/~chenlong/226/Ch6MG.pdf Introduction to Multigrid method> to implement the non-recrusive V-cycle.
%
% _Be careful on the boundary nodes. Restrict smoothing to interiori nodes
% only and enforce the residual on boundary nodes to be zero._
%
% * For one single grid, say J = 4, show the decrease of the residual in
% certain norm for each iteration of the multigrid method.
%
% * Test V-cycle MG for J = 3:6. List iteration steps and cpu time.

%% *Algebraic Multigrid Method*
% 
% We consider solving an SPD matrix equation |Ax = b|, where |A| could be
% obtained as the finite element discretization on a unstructured grids. A
% coarsening of the graph of A is needed and restriction and prolongation
% can be constructued based on the coarsening. A breif introduction on AMG
% can be found at <https://e-reports-ext.llnl.gov/pdf/333205.pdf Falgout, RD. An Introduction to Algebraic Multigrid>

%% Step 1 Matrix
%
% * Load the |lakemesh.mat| in ifem/data.
%
% * Assemble the stiffness matrix on this mesh and take the submatrix
% associated to interior nodes only.
%
% _The mesh is only used to generate the matrix. In the later step, only the
% generated matrix is used._

load lakemesh.mat;
figure; clf; showmesh(node,elem);
A = assemblematrix(node,elem);

%% Step 2 Coarsening
% Use |coarsenAMGc| to get a set of coarse nodes.
help coarsenAMGc

%% Step 3 Transfer Operator and Smoothers
%
% Use |interpolationAMGs| to get restriction and prolongation operators.
help interpolationAMGs

%% Step 4 V-cycle Multigrid
%
% This step is exact the Step 3 in part 3. To test the robustness, apply
% |uniformrefine| to the mesh and generate corresponding matrix.
