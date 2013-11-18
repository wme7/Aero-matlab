%% Equation: Poisson Equation Discretized by $RT_0$ Element in 3D
% We explain the assembling of the matrix equation for the lowest order BDM element
% discretization of Poisson equation. 
%
% [u,sigma] = Poisson3RT0(node,elem,bdEdge,f,g_D,varargin)
%
% Created by Ming Wang at Dec 30. 2010.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Data Structure
%   [elem2dof,dofSign,face] = dof3RT0(elem);
%
% will construct local to global index map; see ifem dof3BDM1doc for
% details.
%
%% Local Bases
% Suppose $i,j,k$ are the vertices of the $l$-th face. The basis are given by
% 
% $$ \phi_l = 2(\lambda_i\nabla \lambda_j \times \nabla\lambda_k +
%             \lambda_j\nabla \lambda_k \times \nabla\lambda_i + 
%             \lambda_k\nabla \lambda_i \times \nabla\lambda_j).$$
%
% Inside one tetrahedron, the 4 bases functions 
% corresponding to 4 local faces [2 3 4; 1 4 3; 1 2 4; 1 3 2] are:
%
% $$ \phi_1 = 2(\lambda_2\nabla \lambda_3 \times \nabla\lambda_4
%             + \lambda_3\nabla \lambda_4 \times \nabla\lambda_2
%             + \lambda_4\nabla \lambda_2 \times \nabla\lambda_3).
% $$
%
% $$ \phi_2 = 2(\lambda_1\nabla \lambda_4 \times \nabla\lambda_3
%             + \lambda_4\nabla \lambda_1 \times \nabla\lambda_3
%             + \lambda_3\nabla \lambda_1 \times \nabla\lambda_4).
% $$
%
% $$ \phi_3 = 2(\lambda_1\nabla \lambda_2 \times \nabla\lambda_4
%             + \lambda_2\nabla \lambda_4 \times \nabla\lambda_1
%             + \lambda_4\nabla \lambda_1 \times \nabla\lambda_2).
% $$
%
% $$ \phi_4 = 2(\lambda_1\nabla \lambda_3 \times \nabla\lambda_2
%             + \lambda_3\nabla \lambda_2 \times \nabla\lambda_1
%             + \lambda_2\nabla \lambda_1 \times \nabla\lambda_3).
% $$
%
% As for the computation of the integration of divergence of the basis,
% We have the following good property.
%
% $$\int_t \phi_i = dofSign(t,i)$$
%
% This property can be drived from relationship between Whitney forms.
% Say the following: 
%
% $$grad w_n = \sum_{e\in patch of n}G_{en}w_e, \qquad
%   curl w_e = \sum_{f\in patch of e}C_{fe}w_f, \qquad
%   div  w_f = \sum_{t\in patch of f}D_{tf}w_t$$
%
% where $w_n, w_e, w_f, w_t$ are nodal element basis, edge element basis, 
% face element basis, and volume element basis respectively.
% As for the volume element (or tetrahedral) basis, 
%
% $$ w_t = 1/v(t) \hbox{ on t and 0 otherwise }.$$
%
% Because of the different oritentation of local and global faces, from
% local bases to the global one, the direction should be corrected. That is
%
% phiGlobal(elem2dof(t,1),:) = phi(t,1)*dofSign(t,1);

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
% jth local faces, suppose i = [i1 i2 i3], j = [j1 j2 j3]. Then
% For two local bases $\phi_i$ and $\phi_j$,
% 
% $$
% \int_T \phi_i \phi_j dx = \int_T \\
% 2(\lambda_{i_1}\nabla \lambda_{i_2}\times\nabla\lambda_{i_3}+
%  \lambda_{i_2}\nabla \lambda_{i_3}\times\nabla\lambda_{i_1}+
%  \lambda_{i_3}\nabla \lambda_{i_1}\times\nabla\lambda_{i_2})\cdot
% $$
% 
% $$
% 2(\lambda_{j_1}\nabla \lambda_{j_2}\times\nabla\lambda_{j_3}+
%  \lambda_{j_2}\nabla \lambda_{j_3}\times\nabla\lambda_{j_1}+
%  \lambda_{j_3}\nabla \lambda_{j_1}\times\nabla\lambda_{j_2})dx
% $$
%
%% Matrix for Differential Operator
% We record $\nabla \cdot \phi_i$ and then the computation $\int _T \nabla
% \cdot \phi_i dx$ is straightforward. Just remember to correct the direction.

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
% and the evulation Jp is just computed once.
%
% The local to global assembling is computed using accumarray
%
%   b = accumarray(elem2dof(:),bt(:),[Ndof 1]);





