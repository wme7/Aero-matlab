
% N = number of points to solve at
N = 500;

% range [a,b) to solve on
a = -1.;
b = 1.;

% discretize the x axis with N points [a,b)
x = linspace(a,b-(b-a)/N, N)';

% compute spacing (here dx is constant)
dx = x(2)-x(1);

% initial condition
u = cos(pi*x);
u = exp(-30*x.*x);

% parameters
U = 2;
nu = 0.1;

% compute time step
dt = 1/( (2*U/dx) + 4*nu/(dx*dx));
dt = 3.2*dt;

% simulation time
StartTime = 0.;
FinalTime = 10.;

% correct dt so that last time step ends at FinalTime
Nsteps = ceil(FinalTime/dt)
dt = FinalTime/Nsteps;

% indices of nodes and left,right offsets
idc = (1:N)';
idl = [N;(1:N-1)'];
idr = [(2:N)';1];

% Jameson-Schmidt-Turkel 
for n=2:Nsteps
    
    v = u;
    for rk=5:-1:1
     % upwind derivative
     dvdx   = (v(idc)-v(idl))/dx;
     d2vdx2 = (v(idl)-2*v(idc)+v(idr))/(dx*dx);

     v = u+(dt/rk)*(-U*dvdx + nu*d2vdx2);
    end
    u = v;
    
    if(~mod(n, 100))
        plot(x,u, 'r'); axis([-1 1 -1 1]); drawnow; pause(0.05);
    end
end

% plot solution at FinalTime
