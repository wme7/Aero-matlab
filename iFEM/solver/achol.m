function [L,p,Ac] = achol(A)

N = size(A,1);
%% Multilvel factorization
level = ceil(log(N)/log(2));
% level = 4;
R = cell(level,1);
isM = cell(level,1);
Nl = zeros(level,1);
Nf = zeros(level,1);
for k = 1:level-1
    Nl(k) = size(A,1);
    [isM{k},As] = md(A);
    [R{k},A] = af(A,As,isM{k});
    Nf(k) = sum(isM{k});
end
Nl(level) = size(A,1);
Nf(level) = Nl(level);
isM{level} = true(Nl(level),1);
R{level} = speye(Nl(level),Nl(level));
% R{level} = spdiags(sqrt(diag(A)),0,Nl(level),Nl(level));
Ac = A;
% R{level} = chol(A);
% opts.type = 'nofill'; opts.michol = 'on';
% R{level} = ichol(A,opts);

%% Ordering/Permutation for the lower part
pl = cell(level,1);
pl{level-1} = 1:Nl(level);
for k = level-2:-1:1
    Ck = (1:Nl(k+1))';     % index for coarse nodes in level k
    Ck1 = Ck(isM{k+1});    % fine nodes in level k+1  goes first
    Ck2 = Ck(~isM{k+1});
    pl{k} = [Ck1; Ck2(pl{k+1})];
           % Ck(isM{k+1})  
           % Ck(P{k+1})    map ordering of level k+1 to level k
    R{k}(:,Ck+Nf(k)) = R{k}(:,pl{k}+Nf(k));
end
coarseNode = find(~isM{1});
p = [find(isM{1}); coarseNode(pl{1})];

%% Merge multilevel R to one lower triangular matrix
L = sparse(N,N);
L(1:N,1:Nf(1)) = transpose(R{1});
col = cumsum(Nf);
for k = 2:level
    L(col(k-1)+1:N,col(k-1)+1:col(k)) = transpose(R{k});
end