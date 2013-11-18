%% Equation: Poisson Equation Discretized by $RT_0$ Element in 2D
% We explain the assembling of the matrix equation for the lowest order 
% Raviart-Thomas element discretization of Poisson equation. 
%
% [u,sigma] = PoissonRT0(node,elem,bdEdge,f,g_D,varargin)
%
% Created by Ming Wang at Jan 2. 2011 and revised to asecnd ordering system
% by Long Chen. Further clearn up by Ming Wang.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[node,elem] = squaremesh([0 1 0 1], 0.5);

%% Data Structure
%
% We use ascend ordering system. Note that the signed area of some
% triangles could be negative.
bdFlag = setboundary(node,elem,'Dirichlet');
[elem,bdFlag] = sortelem(elem,bdFlag);
showmesh(node,elem);
findnode(node);
findelem(node,elem);
display(elem);
%%
% The three local edges are |locEdge = [2 3; 1 3; 1 2]|. The pointer from
% local to global index map can be constructured by
[elem2dof,edge] = dofedge(elem);
findedge(node,edge);
display(elem2dof);
%%
% The global and local orientation of edges are induced from the ascend ordering of
% vertices. Thanks to the ascend ordering, the local and global orientation
% are consistent.
%
% However the ordering orientation is not consistent with the induced
% orientation. More specially the second edge would be [3 1] for consistent
% orientation. So |[1 -1 1]| is used in the construction of div operator.

%% Local bases of RT0 element
%
% Suppose [i,j] is the k-th edge. The two dimensional curl is a rotated
% graident defined as $\nabla^{\perp} f = (-\partial_y f, \partial _x f).$
% The basis of this edge along with its divergence are given by
% 
% $$ \Phi_k = \lambda_i \nabla^{\perp} \lambda_j - \lambda_j \nabla^{\perp} \lambda_i. $$
%
% Inside one triangular, the 3 bases along with their divergence (*where is
% the divergence?*) corresponding to 3 local edges [2 3; 1 3; 1 2] are:
%
% $$ \Phi_1 = \lambda_2 \nabla^{\perp} \lambda_3 - \lambda_3 \nabla^{\perp} \lambda_2. $$ 
%
% $$ \Phi_2 = \lambda_1 \nabla^{\perp} \lambda_3 - \lambda_3 \nabla^{\perp} \lambda_1. $$
%
% $$ \Phi_3 = \lambda_1 \nabla^{\perp} \lambda_2 - \lambda_2 \nabla^{\perp} \lambda_1. $$
%
% The dual basis is the line integral over an orientated edge
%
% $$\int_{e_i} \phi_j de_i = \delta(i,j).$$ 

%% Local bases of P0 element
%
% For triangle t, the basis for the constant function space is p = 1. So in
% the computation of divergence operator, elemSign should be used.

%% Mass Matrix
%
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
% For two local bases $\Phi_i$ and $\Phi_j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2].
%
% $$\int_T \Phi_i \Phi_j dx = \int_T (
%  \lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1}
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1})dx
% $$

%% Matrix for divergence operator
% More here.
[Dlambda,area,elemSign] = gradbasis(node,elem);
B = icdmat(double(elem2dof),elemSign*[1 -1 1]);

%% Right hand side
%
%  The basis for pressure contains a sign. 

%% Boundary condition
%
% use bdFlag to find boundary edges and compute edgeinterpolate with normal
% direction.
