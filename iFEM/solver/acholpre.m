function x = acholpre(b,A,L,R,p,Ac)

N = size(A,1);
Nc = size(Ac,1);

b = b(p);

%% Pre-smoothing by Gauss-Seidel
% x = zeros(N,1);
x = tril(A)\b;
r = b - A*x;

%%
e = L\r;
e(N-Nc+1:N) = Ac\e(N-Nc+1:N); % exact solve in the coarse grid 
e = R\e;      % solve by achol decomposition L*L'
x = x + e;

%% Post-smoothing
x = x + triu(A)\(r-A*e);
x(p) = x;          % permutation back