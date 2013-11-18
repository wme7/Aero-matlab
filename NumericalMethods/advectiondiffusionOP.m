
% N = number of points to solve at
N = 20;

% range [a,b) to solve on
a = -1.;
b = 1.;

% discretize the x axis with N points [a,b)
x = linspace(a,b-(b-a)/N, N)';

% compute spacing (here dx is constant)
dx = x(2)-x(1);

% parameters
U = 2;
nu = 1;

% operator   U*dudx + nu*d2udx2
alpha = U/dx;
beta  = nu/(dx*dx);

d = ones(N,1)*[alpha+beta, -alpha-2*beta, beta];
L = spdiags(d, [-1 0 1], N, N);
L(1,N) = alpha+beta;
L(N,1) = beta;

evs = sort(eig(full(L)))

maxcomputed = max(abs(evs))
maxpredicted =  4*beta+2*alpha