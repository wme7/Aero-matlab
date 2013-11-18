%% Equation: Maxwell Equation Discretized by Edge Element
% We explain the assembling of the matrix equation for the lowest order edge element
% discretization of Maxwell equation. 
%
% [u,edge,A,M,b] = Maxwell(node,elem,bdFace,mu,omega,J)
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Data Structure
%    [elem2dof,dofSign,edge] = dofedge(elem);
%
% will construct local to global index map; see ifem dof3edgedoc for
% details.
%
%% Local Bases
% Suppose [i,j] is the kth edge and i<j. The basis is given by 
% 
% $$ \phi _k = \lambda_i\nabla \lambda_j - \lambda_j \nabla \lambda_i,\qquad
%    \nabla \times \phi_k = 2\nabla \lambda_i \times \nabla \lambda_j.$$
%
% Inside one tetrahedron, the 6 bases functions along with their curl
% corresponding to 6 local edges [1 2; 1 3; 1 4; 2 3; 2 4; 3 4] are
%
% $$ \phi_1 = \lambda_1\nabla\lambda_2 - \lambda_2\nabla\lambda_1,\qquad
%    \nabla \times \phi_1 = 2\nabla\lambda_1\times \nabla\lambda_2,$$
%
% $$ \phi_2 = \lambda_1\nabla\lambda_3 - \lambda_3\nabla\lambda_1,\qquad
%    \nabla \times \phi_2 = 2\nabla\lambda_1\times \nabla\lambda_3,$$
%
% $$ \phi_3 = \lambda_1\nabla\lambda_4 - \lambda_4\nabla\lambda_1,\qquad
%    \nabla \times \phi_3 = 2\nabla\lambda_1\times \nabla\lambda_4,$$
%
% $$ \phi_4 = \lambda_2\nabla\lambda_3 - \lambda_3\nabla\lambda_2,\qquad
%    \nabla \times \phi_4 = 2\nabla\lambda_2\times \nabla\lambda_3,$$
%
% $$ \phi_5 = \lambda_2\nabla\lambda_4 - \lambda_4\nabla\lambda_2,\qquad
%    \nabla \times \phi_5 = 2\nabla\lambda_2\times \nabla\lambda_4,$$
%
% $$ \phi_6 = \lambda_3\nabla\lambda_4 - \lambda_4\nabla\lambda_3,\qquad
%    \nabla \times \phi_6 = 2\nabla\lambda_3\times \nabla\lambda_4.$$
%
% Because of the different oritentation of local and global edges, from
% local bases to the global one, the direction should be corrected. That is
%
%    phiGlobal(elem2dof(t,1),:) = phi(t,1)*dofSign(t,1);

%% Mass Matrix
% We use the integral formula 
%  
% $$ \int_T
% \lambda_1^{\alpha_1}\lambda_2^{\alpha_2}\lambda_3^{\alpha_3}\lambda_4^{\alpha_4}
% dx = \frac{\alpha_1!\alpha_2!\alpha_3!\alpha_4!3!}{(\sum _{i=1}^4\alpha_i
% + 3)!}\;|T|,$$
%
% to get 
% 
% $$ \int _T\lambda_i\lambda_j dx = (1+(i==j))|T|/20. $$
%
% For two local bases $\phi _i$ and $\phi _j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2]. Then
%
% $$ \int _T \phi _i \phi _j dx = \int _T (
% \lambda _{i_1}\lambda _{j_1} \nabla \lambda _{i_2}\nabla\lambda _{j_2} -
% \lambda _{i_1}\lambda _{j_2} \nabla \lambda _{i_2}\nabla\lambda _{j_1} -
% \lambda _{i_2}\lambda _{j_1} \nabla \lambda _{i_1}\nabla\lambda _{j_2} +
% \lambda _{i_2}\lambda _{j_2} \nabla \lambda _{i_1}\nabla\lambda _{j_1} )
% dx. $$

%% Matrix for Differential Operator
% We record $\nabla \times \phi_i$ and then the computation $\int _T \nabla
% \times \phi_i \cdot \nabla \times \phi_j dx$ is straightforward. Just
% remember to correct the direction.

%% Right hand side
% We use 5-points quadrature which is exact for cubic polynomials. In the
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

%%
% Note that the two for loops are nested in such a way that the point pxy
% and the evulation Jp is just computed and stored once.
%
% The local to global assembling is computed using accumarray
%
%   b = accumarray(elem2dof(:),bt(:),[Ndof 1]);

%% Dirichlet boundary condition
% 1. find boundary dof
% 2. assign values to boundary dof




