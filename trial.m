

% Initial Parameters
cfl    = 0.5;
u      = 0.08;
deltax = 1/50; % at least 51 grid points
deltat = cfl*deltax/u;
t_end  = 8 ;

% grid
x = 0:deltax:1;
t = 0:deltat:t_end;
X = length(x);
N = length(t);
T = zeros(N,X);

% Initial Condition
T(1,1:12)  = 1.-(10.*x(1:12)-1).^2;
T(1,12:51) = 0;

% Numerical Scheme:
for  n=1:N-1
    T(n+1,1) = 0; T(n+1,X) = 0;
    for j=2:X-1
    T(n+1,j) = T(n,j) + 1/2*cfl*(T(n,j+1) - T(n,j-1));
    end
end

plot(T(N,:),'red')
plot(T(N,:),'red')