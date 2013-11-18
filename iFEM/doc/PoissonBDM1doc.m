%% Equation: Poisson Equation Discretized by $BDM_1$ Element in 2D
% We explain the assembling of the matrix equation for the lowest order BDM element
% discretization of Poisson equation. 
%
% [u,sigma] = PoissonBDM1(node,elem,bdEdge,f,g_D,varargin)
%
% Created by Ming Wang at Jan 2. 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Data Structure
%   [elem2dof,dofSign,edge] = dofRT0(elem);
%
% will construct local to global index map; see ifem dofBDM1doc for
% details.
%
%% Local Bases
% Suppose [i,j] is the k-th edge. The hierarchic basis along with
% their div are given by
% 
% $$ \phi_k = \lambda_i rot \lambda_j - \lambda_j rot \lambda_i,\quad
%    \psi_k = \lambda_i rot \lambda_j + \lambda_j rot \lambda_i. $$
%
% $$ div \phi_k = 2 \nabla \lambda_i \cdot rot \lambda_j, \quad div \psi_k = 0.$$
%
% Inside one triangular, the 6 bases functions along with their div
% corresponding to 3 local edges [2 3; 3 1; 1 2] are:
%
% $$ \phi_1 = \lambda_2 rot \lambda_3 - \lambda_3 rot \lambda_2,\quad
%    \psi_1 = \lambda_2 rot \lambda_3 + \lambda_3 rot \lambda_2,\quad
%    div\phi_1 = 2 \nabla \lambda_2 \cdot rot \lambda_3, \quad div \psi_1 = 0. $$
%
% $$ \phi_2 = \lambda_3 rot \lambda_1 - \lambda_1 rot \lambda_3,\quad
%    \psi_2 = \lambda_3 rot \lambda_1 + \lambda_1 rot \lambda_3,\quad
%    div\phi_2 = 2 \nabla \lambda_3 \cdot rot \lambda_1, \quad div \psi_2 = 0. $$
%
% $$ \phi_3 = \lambda_1 rot \lambda_2 - \lambda_2 rot \lambda_1,\quad
%    \psi_3 = \lambda_1 rot \lambda_2 + \lambda_2 rot \lambda_1,\quad
%    div\phi_3 = 2 \nabla \lambda_1 \cdot rot \lambda_2, \quad div \psi_3 = 0. $$
%
% Locally, we order the local bases in the following way: 
% 
% $$\{\phi_1,~\,\phi_2,~\,\phi_3,~\,\psi_1,~\,\psi_2,~\, \psi_3.\}$$
%
% Note that $RT_0 \subset BDM_1$, and $\{ \phi_1,~\,\phi_2,~\,\phi_3 \}$ is the 
% local bases for $RT_0$.
%
% Because of the different oritentation of local and global faces, from
% local bases to the global one, the direction should be corrected. That is
%
%    phiGlobal(elem2dof(t,1),:) = phi(t,1)*dofSign(t,1);

%% Mass Matrix
% We use the integral formula 
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
% For two local bases $\phi_i$ and $\phi_j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2].
% $$\int_T \phi_i \phi_j dx = \int_T (
%  \lambda_{i1} rot \lambda_{i2}\cdot\lambda_{j1} rot \lambda_{j2}
% -\lambda_{i2} rot \lambda_{i1}\cdot\lambda_{j1} rot \lambda_{j2}
% -\lambda_{i1} rot \lambda_{i2}\cdot\lambda_{j2} rot \lambda_{j1}
% +\lambda_{i2} rot \lambda_{i1}\cdot\lambda_{j2} rot \lambda_{j1})dx
% $$
%
% $\int_T \psi_i \psi_j dx$ and $\int_T \psi_i \phi_j dx$ can be computed similarly.
%
%% Matrix for Differential Operator
% Note that $\nabla \cdot \psi = 0$, We only need to record $\nabla \cdot \phi_i$ and 
% then the computation $\int _T \nabla
% \cdot \phi_i dx$ is straightforward. Just remember to correct the direction.

%% Right hand side
% We use 5-points quadrature which is exact for cubic polynomials. In the
% barycentric coordinate, the 5-points are
%
% $$ p_1 = [2/3, 1/6, 1/6], \quad  w_1 = 1/3 $$
%
% $$ p_2 = [1/6, 2/3, 1/6], \quad  w_2 = 1/3 $$
%
% $$ p_3 = [1/6, 1/6, 2/3], \quad  w_3 = 1/3 $$
%
%%
% Note that the two for loops are nested in such a way that the point pxy
% and the evulation Jp is just computed once.
%
% The local to global assembling is computed using accumarray
%
%   b = accumarray(elem2dof(:),bt(:),[Ndof 1]);





