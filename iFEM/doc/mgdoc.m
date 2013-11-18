%% MULTIGRID on BISECTION GRIDS
% We describe a geometric-algebraic multigrid methods on bisection grids.
% For example, x = mg(A,b,elem) attempts to solve the system of linear
% equations A*x = b for x. Comparing with algebraic multigrid, we need
% additional information on the linear system: 1. the finest mesh; 2.
% discretization scheme (P1 element). Essentially, it is a geometric
% multigrid. But comparing with geometric multigrid, we do not need the
% hierarchical of the mesh. Instead, we shall call coarsening algorithm to
% recovery the hierarchical structure. This simplifies the data structure
% for adaptive grids and the usage of multigrid methods. Furthermore, the
% hierarchical structure recovered by our coarsening algorithm leads to a
% more efficient smoother and better convergent rate of multigrid or its
% variants.
%
% We should mention that it only works for bisection grids stored in a
% suitable order. See iFEMdoc algorithm/coarsen for details.

%%
% Load A, b, and elem.
load MGdoc
N = length(b); x0 = zeros(N,1);

%% Hierarchical Structure of Mesh
% We use uniformcoarsen to reconstruct a hieratchical structure (HB) of the
% mesh (elem). In HB matrix, HB(:,2:3) are two parent nodes of the node
% HB(:,1), i.e.,
%
%               HB(:,2) --- HB(:,1) --- HB(:,3)
%
% HB(:,1) is the index in the fine level and HB(:,2:3) are in the coarse
% level. The array NL records the range of indices in each level. The nodal
% rangle in the k-th level is given by NL(k)+1:NL(k+1).

HB = zeros(N,3);
level = max(min(round(log2(N))-8,12),2);
NL(level+1) = N; 
NL(1) = 0; % at least two level methods
for k = level: -1 : 2
    [elem,newHB] = uniformcoarsen(elem);
    if isempty(newHB)      % no nodes are removed in the coarsening step
        NL = NL(k:end);       
        break; 
    end
    NL(k) = NL(k+1) - size(newHB,1);
    HB(NL(k)+1:NL(k+1),1:3) = newHB(:,1:3);
end
display(NL);
%%
% One nice feature of our coarsening algorithm is the number of nodes is
% almost reduced by half since nodes added in different refinement step
% could be coarsened at the same time.

%% Transfer Operators
% Construct transfer operators:
%%
% 
% * Pro: prolongation operator from coarse space to fine space
% * Res: restriction operator from fine space to coarse space
% * Ai: matrix in each level 
% * Bi: iterator(smoother) in each level

level = length(NL)-1; 
Pro = cell(level,1);          
Res = cell(level,1);          
Ai = cell(level,1);           
Bi = cell(level,1);           
BBi = cell(level,1);          
Ai{level} = A;                
for j = level:-1:2
    fineNodeRange = NL(j)+1:NL(j+1);
    fineNode = HB(fineNodeRange,1);
    coarseNode = (1:NL(j))';
    nFineNode = NL(j+1)-NL(j);
    isCoarseNode = true(NL(j+1),1);
    isCoarseNode(fineNode) = false;
    coarseNodeFineIdx = find(isCoarseNode);
    ii = [coarseNodeFineIdx; fineNode; fineNode];
    jj = [coarseNode; HB(fineNodeRange,2); HB(fineNodeRange,3)];
    ss = [ones(NL(j),1); 0.5*ones(nFineNode,1); 0.5*ones(nFineNode,1)];
    Pro{j-1} = sparse(ii,jj,ss,NL(j+1),NL(j));
    Res{j} = Pro{j-1}';         % restriction is the transpose of prologation              
    Ai{j-1} = Res{j}*Ai{j}*Pro{j-1};  % Ac = Res*Af*Pro
    Bi{j} = tril(Ai{j});        % presmoothing B = D+L
    BBi{j} = triu(Ai{j});       % postsmoothing BB = D+U
end
%%
[Prolongation,Restriction] = prolongationexample;
%% 
% Left: fine mesh; Right: coarse mesh.
%%
% * fine node: 1,2,3,4,5,6,7
% * coarse node: 1,2,3,4,5
%
% Note that the 5-th node in the coarse mesh is the 7-th node in the fine
% mesh. The array coarseNodeFineIdx is to find the indices of coarse nodes
% in the fine mesh.
%
% For this simple example
%
% * NL = [0 5 7]; 
% * HB(6,1) = 5; HB(6,2) = 1; HB(6,3) = 4;
% * HB(7,1) = 6; HB(7,2) = 3; HB(7,3) = 4;
% * fineNode = [5 6];
% * coarseNode = [1 2 3 4 5];
% * coarseNodeFineIdx = [1 2 3 4 7];

%%
makeHtmlTable(Prolongation);
%%
% Prolongation matrix
makeHtmlTable(Restriction);
%%
% Restriction matrix

%% V-cycle Multigrid
% One iteration x = x + B*(f-A*x) consists of three steps:
%%
% 
% # Form residual: r = b-A*x;
% # Compute approximate correction: e = B*r;
% # Correct the solution: x = x + e;
%
% The action B*r is computed using V-cycle multigrid.

r = cell(level,1);            % residual in each level
e = cell(level,1);            % correction in each level
k = 1;
tol = 1e-12;
x = x0;
ek = ones(N,1);
error = norm(b-A*x)/norm(b);
while error > tol
    k = k+1;
    % Step 1: Form residual r
    rk = b - A*x;       
    % Step 2: Compute B*r by MG
    r{level} = rk;
    for i = level:-1:2
        e{i} = Bi{i}\r{i};  % pres-moothing
        r{i-1} = Res{i}*(r{i}-Ai{i}*e{i}); % restriction of the residual
    end
    e{1} = Ai{1}\r{1};      % exact solver in the coarsest mesh
    for i = 2:level
        e{i} = e{i} + Pro{i-1}*e{i-1};     % prolongation of the correction
        e{i} = e{i} + BBi{i}\(r{i}-Ai{i}*e{i}); % post-smoothing
    end
    ekold = ek;
    ek = e{level};
    % Step 3: Correct the solution
    x = x + ek;
    error = norm(ek)/norm(x);
%    error = norm(b-A*x)/norm(b);
end
fprintf('Number of unknowns: %8.0u, MG iteration: %2.0u, error = %12.8g\n',...
         length(b),k,error)
     
%% MultiGrid Preconditioner
% The V-cycle operator B is symmetric and positive definite. We can use it
% as a preconditioner and PCG (preconditioned Conjugate Gradient) method.
% 
% Given $x_0; r_0 = b - Ax_0; p_0 = Br_0$
%
% for $k=1,...,$ till convergence
% 
% $$\alpha_k = (Br_{k-1},r_{k-1})/(Ap_{k-1},p_{k-1});$$
% 
% $$r_k = r_{k-1} - \alpha _k A p_{k-1};$$
%
% $$x_k = x_{k-1} + \alpha _kp_{k-1};$$
%
% $$\beta_k = (Br_{k},r_{k})/(Br_{k-1},r_{k-1});$$
%
% $$p_k = Br_k + \beta_k p_{k-1};$$
%
% end
r = cell(level,1);            % residual in each level
e = cell(level,1);            % error in each level
k = 1;
x = x0;
rk = b - A*x;
error = norm(b-A*x)/norm(b);
while error > tol
    % compute Br by MG
    r{level} = rk;
    for i=level:-1:2
        e{i} = Bi{i}\r{i};
        r{i-1} = Res{i}*(r{i}-Ai{i}*e{i});
    end
    e{1} = Ai{1}\r{1};         % direct solver in the coarest level  
    for i=2:level
        e{i} = e{i} + Pro{i-1}*e{i-1};
        e{i} = e{i} + BBi{i}\(r{i}-Ai{i}*e{i});
    end
    Brk = e{level};
    % update tau, beta, and p
    tauk = Brk'*rk;
    if k==1
        pk = Brk;
    else
        beta = tauk/tauold;
        pk = Brk + beta*pk;
    end
    tauold = tauk;
    % update alpha, x, and r
    Apk = A*pk;
    alpha = tauk/(Apk'*pk);
    rk = rk - alpha*Apk;
    x = x + alpha*pk;
    % compute error for the stopping criterion
%     error = norm(alpha*pk)/norm(x);  % relative error of two iterations 
    error = alpha*sqrt(pk'*Apk/(x'*A*x)); % relative error in energy norm
%     error = norm(b-A*x)/norm(b);  % relative error of the residual
    k = k + 1;
end
fprintf('Number of unknowns: %8.0u, MGCG iteration: %2.0u, error = %12.8g\n',...
         length(b),k, error)
%%
% MGCG (multigrid as a preconditioner in CG) converges faster than
% multigrid itself. It is also more robust. For example, MGCG converges
% uniformly with respect to the jump of diffusion coefficients while MG
% alone may not be.

%% Remarks
% 1. For 3D bisection grids, the only difference is the reconstruction of
% hierarchical structure using uniformcoarsen3. See MG3P1 and
% uniformcoarsen3 for details.
%
% 2. For P2 elements, we simply add one level in the construction of
% hierarchical structure. Namely nodes in the finest level are degree of
% freedom associated to edges. See MGP2, MG3P2 for details.
%
% 3. For adaptive grids, to acheive optimal computational cost, the
% smoothing should be restricted to the new added nodes and two parent nodes. 
% It is therefore called local multigrid. Since our coarsening algorithm
% will produce a sequence of nested meshes whose number of vertices is
% geometrically decreased by a factor 1/2, we do not implement local MG.



