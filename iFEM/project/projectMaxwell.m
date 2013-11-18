%% Project: Edge Finite Element Method for Maxwell-type Equations
%
% The purpose of this project is to implement the lowest order edge element
% for solving the following eddy current equations in two dimensions
%
% $$ \nabla \times \nabla \times \mathbf u +  \mathbf u = \mathbf f $$
%
% with Dirichlet boundary condition $$ \mathbf u\times \mathbf n = g_D. $$

%% Step 1: Data Structure
%
% Since the unknowns are associated to edges, you need to generate edges
% and more importantly the indices map from a triangle to global indices of
% its three edges.  In each triangle, the three edges are labled such that
% the i-th edge is opposite to the i-th vertex for i=1,2,3 in each triangle
% and oritentated counterclockwise i.e., [2 3; 3 1; 1 2]. The orientation
% of the edge in the triangulation is from the node with smaller global
% index to bigger one. The output dofSign records the consistency of the
% local and global edge orientation.

[elem2edge,edge,dofSign] = dofedge(elem);

%% Step 2: Assemble matrices
%
% Suppose [i,j] is the kth edge and orientated from i to j. The basis is
% given by
% 
% $$ \phi _k = \lambda_i\nabla \lambda_j - \lambda_j \nabla \lambda_i.$$
%
% Use the following subroutine to generate $$ \nabla \lambda_i $$

[Dlambda,area] = gradbasis(node,elem);
%%
% and compute $$ curl \phi _k $$ which is a piecewise constant. Then the entry
% $$ ( curl \phi_k, curl \phi_l) $$ can be computed accordingly. When
% assemble the local entry to global entry, don't forgot the sign
% consistency.

%% 
% Use the integral formula 
%  
% $$ \int_T
% \lambda_1^{\alpha_1}\lambda_2^{\alpha_2}\lambda_3^{\alpha_3}
% dx = \frac{\alpha_1!\alpha_2!\alpha_3!2!}{(\sum _{i=1}^3\alpha_i
% + 2)!}\;|T|,$$
%
% to get 
% 
% $$ \int _T\lambda_i\lambda_j dx = (1+(i==j))|T|/12. $$
%
% Use the above formula to compute the mass matrix $$ (\phi_k, \phi_l). $$
% Again don't forget to correct the sign.

%% Step 3: Right hand side
%
% Use one point quadrature to compute $$ (f, \phi_k). $$ Remember to
% correct the sign.

%% Step 4: Modify equations to include boundary conditions
%
% Modify the right hand side to include Dirichlet boundary
% conditions. 
%
% To find out the boundary edges, use

s = accumarray(elem2edge(:), 1, [size(edge,1) 1]);

%%
% The boundary value associated to edges on the boundary is given by the
% edge integral
%
% $$ \frac{1}{|E|}\int_E \mathbf u \cdot \mathbf t \, ds. $$
%
% Use Simpson rule to compute an approximation of the above integral.

%% Step 5: Solve the linear system
%
% Use direct solver to solve the linear system.

%% Step 6: Verify the convergence 
%
% - Substitude a smooth function into the equation to get a right hand side.
% Pass the data f and g_D to your code to compute an approximation u.
%
% - Compute the edge interpolant u_I by computing edge integrals using the
% exact formula of u.
%
% - Compare u_I and u in the energy norm given by the problem.
%
% - Refine mesh several times to show the order of convergences.
