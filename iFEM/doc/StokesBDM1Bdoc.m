%% Equation: Stokes Equations Discretized by $RT_0$ Element in 2D
% We explain the assembling of the matrix equation for the lowest order 
% BDM element with bubble discretization of Stokes equations. 
%
% [u,p,w,rotuI,Me,Mv,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option)
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
elemunSort = elem;
[elem,bdFlag] = sortelem(elemunSort,bdFlag);
showmesh(node,elem);
findnode(node);
findelem(node,elem);
display(elem);
%%
% The three local edges are |locEdge = [2 3; 1 3; 1 2]|. The pointer from
% local to global index map can be constructured by
[elem2edge,edge] = dofedge(elem);
edge = double(edge); elem = double(elem);
elem2edge = double(elem2edge);
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
[Rlambda,signedarea,elemSign] = curlbasis(node,elem);
area = abs(signedarea);
%% Local bases of P2 element with bubble
% The bubble function is:
% $$ w_b  = 27 \lambda_1 \lambda_2 \lambda_3$$
% Inside one triangular, the 6 bases corresponding to 3 local edges [2 3; 1 3; 1 2] are:
% $$ \Phi_1 = (2\lambda_1 - 1)\lambda_1 +1/9 w_b.$$ 
%
% $$ \Phi_2 = (2\lambda_2 - 1)\lambda_2 +1/9 w_b.$$ 
%
% $$ \Phi_3 = (2\lambda_3 - 1)\lambda_3 +1/9 w_b.$$ 
%
% $$ \Phi_4 = 4\lambda_2\lambda_3 - 4/9 w_b. $$ 
%
% $$ \Phi_5 = 4\lambda_1\lambda_3 - 4/9 w_b. $$
%
% $$ \Phi_6 =4\lambda_1 \lambda_2 - 4/9 w_b. $$
%
% The interior bubble basis is
% $$ \Phi_7 = w_b.$$
%
%% Local bases of BDM1 element with bubble
%
% Suppose [i,j] is the k-th edge. The two dimensional curl is a rotated
% graident defined as $\nabla^{\perp} f = (-\partial_y f, \partial _x f).$
% The basis of this edge along with its divergence are given by
% 
% $$ \Phi_k = \lambda_i \nabla^{\perp} \lambda_j - \lambda_j \nabla^{\perp} \lambda_i. $$
%
% Inside one triangular, the 6 bases corresponding to 3 local edges [2 3; 1 3; 1 2] are:
% $$ \Phi_1 = \lambda_2 \nabla^{\perp} \lambda_3 - \lambda_3 \nabla^{\perp} \lambda_2. $$ 
%
% $$ \Phi_2 = \lambda_1 \nabla^{\perp} \lambda_3 - \lambda_3 \nabla^{\perp} \lambda_1. $$
%
% $$ \Phi_3 = \lambda_1 \nabla^{\perp} \lambda_2 - \lambda_2 \nabla^{\perp} \lambda_1. $$
%
% $$ \Phi_4 = \lambda_2 \nabla^{\perp} \lambda_3 + \lambda_3 \nabla^{\perp} \lambda_2. $$ 
%
% $$ \Phi_5 = \lambda_1 \nabla^{\perp} \lambda_3 + \lambda_3 \nabla^{\perp} \lambda_1. $$
%
% $$ \Phi_6 = \lambda_1 \nabla^{\perp} \lambda_2 + \lambda_2 \nabla^{\perp} \lambda_1. $$
%
% The interior bubble basis is
% $$ \Phi_7 = \nabla^{\perp} w_b.$$

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
% The lumped mass matrix of the $P_2$ with bubble function is a diagonal
% matrix. Inside one triangle, the lumped matrix is:
% $$ \int _T\lambda_i\lambda_i dx =|T|/20, i  = 1:3 $$
% $$ \int _T\lambda_i\lambda_i dx =|T|2/15, i = 4:6 $$
% $$ \int _T\lambda_i\lambda_i dx =|T|9/20, i = 7 $$

N = size(node,1); NT = size(elem,1); NE = size(edge,1); 
Nu = 2*NE+NT; Np = NT; 

[Mvb,invMv] = asmMassMatP2b(node,elem);
Mv = asmMassMatP2(node,elem,'NB');

%%
%
% For two local bases $\Phi_i$ and $\Phi_j$ corresponding to the ith and
% jth local edges, suppose i = [i1 i2], j = [j1 j2].
%
% The mass matrix of the BDM1B element:
% 
% For i = 1:3, j = 1:3, 
%
% $$\int_T \Phi_i \Phi_j dx = \int_T (
%  \lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1}
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1})dx
% $$
%
% $$=|T|/12(
%  (1+(i1==j1))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j2}
% -(1+(i2==j1)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j2}
% -(1+(i1==j2))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j1}
% +(1+(i2==j2)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j1})
% $$
%
% For i = 4:6, j = 4:6,
%
% $$\int_T \Phi_i \Phi_j dx = \int_T (
%  \lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% +\lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1}
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1})dx
% $$
%
% $$=|T|/12(
%  (1+(i1==j1))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j2}
% +(1+(i2==j1)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j2}
% +(1+(i1==j2))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j1}
% +(1+(i2==j2)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j1})
% $$
%
% For i = 1:3, j = 4:6 or j = 1:3, i = 4:6
%
% $$\int_T \Phi_i \Phi_j dx = \int_T (
%  \lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j1} \nabla^{\perp} \lambda_{j2}
% +\lambda_{i1} \nabla^{\perp} \lambda_{i2}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1}
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\lambda_{j2} \nabla^{\perp} \lambda_{j1})dx
% $$
%
% $$=|T|/12(
%  (1+(i1==j1))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j2}
% -(1+(i2==j1)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j2}
% +(1+(i1==j2))\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{j1}
% -(1+(i2==j2)) \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{j1})
% $$
%
% For i = 1:3, j = 7 or j = 7, i = 1:3
%
% $$\int_T \Phi_i \Phi_7 dx = \int_T (
%  \lambda_{i1} \nabla^{\perp}\lambda_{i2}\cdot \nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})
% -\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})
% $$
%
% $$=|T|9/10(
% \nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{i3}
% + \nabla^{\perp} \lambda_{i2}\cdot\nabla^{\perp} \lambda_{i2}
% -\nabla^{\perp} \lambda_{i1}\cdot \nabla^{\perp} \lambda_{i3}
% -\nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{i1})
% $$
% 
% For i = 4:6, j = 7 or j = 7, i = 4:6
%
% $$\int_T \Phi_i \Phi_7 dx = \int_T (
%  \lambda_{i1} \nabla^{\perp}\lambda_{i2}\cdot \nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})
% +\lambda_{i2} \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})
% $$
%
% $$=|T|9/10(
% \nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{i3}
% + \nabla^{\perp} \lambda_{i2}\cdot\nabla^{\perp} \lambda_{i2}
% +\nabla^{\perp} \lambda_{i1}\cdot \nabla^{\perp} \lambda_{i3}
% +\nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{i1})
% $$
% 
% For i = 7, j = 7 
%
% $$\int_T \Phi_7 \Phi_7 dx = \int_T (
%  \nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})\cdot \nabla^{\perp}(\lambda_{i1}\lambda_{i2}\lambda_{i3})
% $$
%
% $$=|T|81/10(
% \nabla^{\perp} \lambda_{i1}\cdot \nabla^{\perp} \lambda_{i1}
% + \nabla^{\perp} \lambda_{i1}\cdot\nabla^{\perp} \lambda_{i2}
% +\nabla^{\perp} \lambda_{i1}\cdot \nabla^{\perp} \lambda_{i3}
% +\nabla^{\perp} \lambda_{i2}\cdot\nabla^{\perp} \lambda_{i2})
% +\nabla^{\perp} \lambda_{i2}\cdot \nabla^{\perp} \lambda_{i3}
% +\nabla^{\perp} \lambda_{i3}\cdot\nabla^{\perp} \lambda_{i3})
% $$
%
Me = asmMassMatBDM1(node,elem);
Meb = asmMassMatBDM1B(node,elem);
%%
%
% The mass matrix of the $P_0$ element:
% $$\int_T 1_{T_i}*1_{T_j} dx = \delta(i,j)|T|.$$
%
% The inverse of the mass matrix is:
%
invMt = spdiags(1./area,0,NT,NT);
%% Matrix for  -divergence operator
RT0B = -icdmat(double(elem2edge),elemSign*[1 -1 1]);
B = [RT0B sparse(NT,NE+NT)];

%% Matrix for curl operator
% Firstly, we assemble the curl operator with the Hieratical bases for P_2.
%
% With this basis, the curl operator is a block diagonal matrix. 
%
% The first NE by N block is the curl operator with P_1-RT0 elements. 
% 
% The second NE by NE block is 4 times identity matrix.
%
% The last NT by NT block is identity matrix.

RT0C = icdmat(double(edge),[-1 1]); 
HC = blkdiag(RT0C,4*speye(NE,NE),speye(NT,NT)); 
% Secondly, we find the transfer operator between the P2b bases and the
% Hieratical bases.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the transform matrix between bases
% B1                                                                     B2
% 1) (2\lambda_i -1)\lambda_i +1/9w_b                 1) \lambda_i
% 2) 4\lambda_i \lambda_j -4/9w_b         <===>  2) 4\lambda_i\lambda_j
% 3) 27\lambda_1\lambda_2\lambda_3                  3) 27\lambda_1\lambda_2\lambda_3  
% B1 = B2*T, where T is the tranfer matrix
% T = | 1      0    0 |
%       |1/2    1    0 |
%       |1/9  4/9   1 |
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2NE = sparse([edge(:,1);edge(:,2)],repmat(1:NE,1,2),1,N,NE);
N2NT = sparse([elem(:,1);elem(:,2);elem(:,3)], repmat(1:NT,1,3),1,N,NT);
NE2NT = sparse([elem2edge(:,1);elem2edge(:,2);elem2edge(:,3)],repmat(1:NT,1,3),1,NE,NT);
T = [speye(N,N) sparse(N,NE) sparse(N,NT);
       -1/2*N2NE'  speye(NE,NE)   sparse(NE,NT);
        1/9*N2NT' -4/9*NE2NT' speye(NT,NT)];
    
%
% At last, the curl operator with P2b bases is:
C = HC*T;

%% Matrix for rot operator
R = invMv*C'*Meb;

%% Matrix for vector Laplacian
A = B'*invMt*B + R'*Mvb*R;

%% Right hand side
%
% To assemble the right hand side, we compute the inner product of the
% right hand side function f and the test function with quadrature. The
% default order is 4. The right hand hand side of the pressure equationas
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
% For the velocity right hand side: fu = fu-A*u-Meb*C*invMv*ubd;
%
% For the pressure right hand side:fg = 0-B*u.