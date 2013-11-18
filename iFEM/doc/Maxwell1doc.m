%% Equation: Maxwell Equation Discretized by Edge Element
% We explain the assembling of the matrix equation for the lowest order
% second family edge element discretization of Maxwell equation.
%
% [u,edge,A,M,b] = Maxwell(node,elem,bdFace,mu,omega,J)
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Similarity and Difference
% The second family share the same data structure with the first one. The
% assembeling of the matrix and right hand side is thus similar. In
% particular, the curl-curl matrix is the same (not exactly identical since
% now the dimension is 2*Ndof by 2*Ndof. The first family is only Ndof by
% Ndof. The mass matrix is extended. Also one needs additional d.o.f when
% assign the Dirichlet boundary condition.

%% Data Structure
%    [elem2dof,dofSign,edge] = edgedofstructure(elem);
%
% will construct local to global index map; see ifem edgedofdoc for
% details.

%% Local Bases
% Suppose [i,j] is the kth edge. The basis for the first family is given by 
% 
% $$ \Phi _k = \lambda_i\nabla \lambda_j - \lambda_j \nabla \lambda_i,\qquad
%    \nabla \times \Phi_k = 2\nabla \lambda_i \times \nabla \lambda_j.$$
%
% The additional 6 bases for the second family are:
%
% $$ \Psi _k = \lambda_i\nabla \lambda_j + \lambda_j \nabla \lambda_i,\qquad
%    \nabla \times \Psi_k = 0.$$
%
% Inside one tetrahedron, the 6 bases functions along with their curl
% corresponding to 6 local edges [1 2; 1 3; 1 4; 2 3; 2 4; 3 4] are
%
% $$ \Psi_1 = \lambda_1\nabla\lambda_2 + \lambda_2\nabla\lambda_1,$$
%
% $$ \Psi_2 = \lambda_1\nabla\lambda_3 + \lambda_3\nabla\lambda_1,$$
%
% $$ \Psi_3 = \lambda_1\nabla\lambda_4 + \lambda_4\nabla\lambda_1,$$
%
% $$ \Psi_4 = \lambda_2\nabla\lambda_3 + \lambda_3\nabla\lambda_2,$$
%
% $$ \Psi_5 = \lambda_2\nabla\lambda_4 + \lambda_4\nabla\lambda_2,$$
%
% $$ \Psi_6 = \lambda_3\nabla\lambda_4 + \lambda_4\nabla\lambda_3.$$

%% Matrix for Differential Operator
% The curl-curl matrix is the same as the first family since $\nabla \times
% \Psi = 0$.

%% Mass Matrix
% We use the integral formula 
%  
% $$ \int _T\lambda_i\lambda_j dx = (1+(i==j))|T|/20. $$
%
% For two local bases $\Psi _i$ and $\Psi _j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2]. Then
%
% $$ \Phi _i = \lambda_{i_1}\nabla \lambda_{i_2} - \lambda_{i_2} \nabla
% \lambda_{i_1},\qquad  \Psi _j = \lambda_{j_1}\nabla \lambda_{j_2} +
% \lambda_{j_2} \nabla \lambda_{j_1}.$$
%
% $$ \int _T \Phi _i \Psi _j dx = \int _T (
% \lambda _{i_1}\lambda _{j_1} \nabla \lambda _{i_2}\nabla\lambda _{j_2} +
% \lambda _{i_1}\lambda _{j_2} \nabla \lambda _{i_2}\nabla\lambda _{j_1} -
% \lambda _{i_2}\lambda _{j_1} \nabla \lambda _{i_1}\nabla\lambda _{j_2} -
% \lambda _{i_2}\lambda _{j_2} \nabla \lambda _{i_1}\nabla\lambda _{j_1} )
% dx. $$
%
% $$ \int _T \Psi _i \Psi _j dx = \int _T (
% \lambda _{i_1}\lambda _{j_1} \nabla \lambda _{i_2}\nabla\lambda _{j_2} +
% \lambda _{i_1}\lambda _{j_2} \nabla \lambda _{i_2}\nabla\lambda _{j_1} +
% \lambda _{i_2}\lambda _{j_1} \nabla \lambda _{i_1}\nabla\lambda _{j_2} +
% \lambda _{i_2}\lambda _{j_2} \nabla \lambda _{i_1}\nabla\lambda _{j_1} )
% dx. $$

%% Right hand side
% We still use 5-points quadrature which is exact for cubic polynomials. In the
% barycentric coordinate, the 5-points are
%
% $$ p_1 = [1/4, 1/4, 1/4, 1/4], \quad  w_1 = -4/5 $$
%
% $$ p_2 = [1/2, 1/6, 1/6, 1/6], \quad  w_2 = 9/20 $$
%
% $$ p_3 = [1/6, 1/2, 1/6, 1/6], \quad  w_3 = 9/20 $$
%
% $$ p_4 = [1/6, 1/6, 1/2, 1/6], \quad  w_4 = 9/20 $$
%
% $$ p_5 = [1/6, 1/6, 1/6, 1/2], \quad  w_5 = 9/20 $$

%% Dirichlet boundary condition



