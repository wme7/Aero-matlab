%% Equation: Stokes Equations Discretized by $RT_0$ Element in 2D
% We explain the assembling of the matrix equation for the lowest order 
% Raviart-Thomas element discretization of Stokes equations. 
%
% [u,p,w,rotuI,Me,Mv,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option)
%
%  Basically, we start from the following equation:
%  ---------------------------------------------------------------
%    - grad div u + curl rot u + grad p  = f   in \Omega   
%                              - div u   = 0   in \Omega   
%                                    u   = g   on \Gamma   
%  ---------------------------------------------------------------
% The bilinear form is:
%  ---------------------------------------------------------------------
% (rot u,rot v) +(div u, div v)  - (p, div v)  = (f,v)  v \in $H_0(div)$ 
%                                 (- div u,q)   = 0      q\in $L_0^2$
%                                           u   = g      on \Gamma   
%  ---------------------------------------------------------------------
% In the matrix form:
%  ---------------------------------------------------------------------
% |R'*Mv*R +B'*invMt*B B'||u| = fu   
% |    B               0 ||p| = 0     
%                            u = g      on \Gamma   
%  ---------------------------------------------------------------------
% Created by Lin Zhong at April. 2013.
% See PoissonRT0doc.m 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 
%
clear all; close all;
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
[elem2edge,edge] = dofedge(elem);
findedge(node,edge);
display(elem2edge);
%%
%
% The global and local orientation of edges are induced from the ascend ordering of
% vertices. Thanks to the ascend ordering, the local and global orientation
% are consistent.
%
% However the ordering orientation is not consistent with the induced
% orientation. More specially the second edge would be [3 1] for consistent
% orientation. So |[1 -1 1]| is used in the construction of div operator.
%
% The ordering will introduce signs to the elem. 
% If the assending order of the elem is consist with the right hand rule,
% the sign of the elem is 1, otherwise -1.
[Dlambda,area,elemSign] = gradbasis(node,elem);
%% Local bases of RT0 element
%
% Suppose [i,j] is the k-th edge. The two dimensional curl is a rotated
% graident defined as $\nabla^{\perp} f = (-\partial_y f, \partial _x f).$
% The basis of this edge along with its divergence are given by
% 
% $$ \Phi_k = \lambda_i \nabla^{\perp} \lambda_j - \lambda_j \nabla^{\perp} \lambda_i. $$
%
% Inside one triangular, the 3 bases corresponding to 3 local edges [2 3; 1 3; 1 2] are:
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
% We use the integral formula 
%  
% $$ \int_T
% \lambda_1^{\alpha_1}\lambda_2^{\alpha_2}\lambda_3^{\alpha_3}
% dx = \frac{\alpha_1!\alpha_2!\alpha_3!2!}{(\sum _{i=1}^3\alpha_i
% + 2)!}\;|T|.$$
%
% The mass matrix of the $P_1$ element:
%
% $$ \int _T\lambda_i\lambda_j dx = (1+(i==j))|T|/12. $$
%  
% We apply mass lumping, and get:
%  $$ \int _T\lambda_i\lambda_i dx =(i==j)|T|/3. $$
%

N = size(node,1); NT = size(elem,1); NE = size(edge,1); 
Nu = NE; Np = NT;
Mv = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
invMv = spdiags(1./Mv,0,N,N);
Mv = spdiags(Mv,0,N,N);

%%
%
% For two local bases $\Phi_i$ and $\Phi_j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2].
%
% The mass matrix of the RT0 element:
%
% $$\int_T \Phi_i \Phi_j dx = \int_T (
%  \lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1}
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1})dx
% $$
%
Me = getmassmatvec(elem2edge,area,Dlambda,'RT0');
%%
%
% The mass matrix of the $P_0$ element:
% $$\int_T 1_{T_i}*1_{T_j} dx = \delta(i,j)|T|.$$
%
% The inverse of the mass matrix is:
%
invMt = spdiags(1./area,0,NT,NT);
%% Matrix for  -divergence operator
B = icdmat(double(elem2edge),elemSign*[1 -1 1]);

%% Matrix for curl operator
C = icdmat(double(edge),[-1 1]);

%% Matrix for rot operator
R = invMv*C'*Me;

%% Matrix for vector Laplacian
A = B'*invMt*B + R'*Mv*R;

%% Right hand side
%
% To assemble the right hand side, we compute the inner product of the
% right hand side function f and the test function with quadrature. The
% default order is 3. The right hand hand side of the pressure equationas
% are zeros.

%% Boundary condition
%
%% 
%
% (1) Use the bdFlag to find the Dirichlet boundary and find the edgeSign of
% the boundary edges.
% About how to construct edgeSign, please see:...
%
% (2) According to the integration by parts formula 
% $$(curl u, tau) = (u, curl \tau) +(u \cdot t, tau)|_{\partial \Omega},$$
% we need to compute the line integral of $u\cdot t$ on the boundary, i.e.,
% ubd = (u\cdot t, \tau)|_{\partial \Omega}.
%
% (3) Modify the matrix and the right hand side. 
%%
% The boundary part of the matrix should be identity matrix.
%%
% 
% For the velocity right hand side: fu = fu-A*u-Me*C*invMv*ubd;
%
% For the pressure right hand side:fg = 0-B*u.