
function fvsolver(N, method) 
% N = number of grid points

% grid spacing on interval (-1,1)
dx = 2./N; 

% location of cell centers 
x = linspace(-1+0.5*dx, 1-0.5*dx, N);

% advection speed
u = 1;

% cfl condition + safety margin
dt = 0.2*dx;

% final time 
FinalTime = 2;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% initial condition
q = fvexact(x); % q = 1*(abs(x)<=0.25);

% time step
for n=1:Nsteps
  
  qp1  = [q(2:N),q(1)];    % q_{m+1}
  qm1 = [q(N),q(1:N-1)]; % q_{m-1}

  % compute slopes s_m based on limiter choice
  
  if(strcmp(method,'minmod'))
      s = minmod((q-qm1)/dx, (qp1-q)/dx);
  elseif(strcmp(method,'MC'))
      s = minmod( (qp1-qm1)/(2*dx), minmod(2*(q-qm1)/dx, 2*(qp1-q)/dx));
  elseif(strcmp(method,'Fromm'))
      s = (qp1-qm1)/(2*dx);
  elseif(strcmp(method,'BeamWarming'))
      s = (q-qm1)/(dx); % upwind
  elseif(strcmp(method,'LaxWendroff'))
      s = (qp1-q)/(dx); % downwind
  else
      s = zeros(size(q));
  end
 
  sm1 = [s(N),s(1:N-1)];                           % s_{m-1}
  
  % update cell averages
  q = q - (u*dt/dx)*(q-qm1) - 0.5*(u*dt/dx)*(dx-u*dt)*(s-sm1); 

  if(mod(n,40)==0)
    plot(x,q, '.'); pause(0.02);
  end
end

plot(x, q, 'b.'); 
% have gone one exact period
hold on; plot(x, fvexact(x), 'r-'); hold off;
xlabel('x'); 
legend(sprintf('%s N=%d T=%f', method, N, FinalTime), 'Exact solution')
axis([-1 1 -.5 1.5])