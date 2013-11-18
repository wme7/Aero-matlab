
function fluxlimiter(N, method) 
% N = number of grid points

% grid spacing on interval (-1,1)
dx = 2./N; 

% location of cell centers 
x = linspace(-1+0.5*dx, 1-0.5*dx, N);

% advection speed
u = 1;

% cfl condition + safety margin
dt = 0.25*dx;

% final time 
FinalTime = 2;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
lam = u*dt/dx;

% initial condition
q = fluxlimiterexact(x); 

% time step
for n=1:Nsteps
  
  qp1  = [q(2:N),q(1)];    % q_{m+1}
  qm1 = [q(N),q(1:N-1)]; % q_{m-1}
  qm2 = [qm1(N), qm1(1:N-1)]; % q_{m-2}
  
  % compute limited cell average jumps based on method choice
  % theta_{i-1/2} = (q_{i-1}-q_{i-2})/(q_{i}-q_{i-1})
  tol = 1e-10;
  dq1 = q-qm1;
  dq2 = qm1-qm2;
  % have to test in case dq1 is small because of finite precision effects
  
  % are both dq1, and dq2 small (smooth)
  ids = find(abs(dq1)+abs(dq2)<tol); % both jumps small
  thetaL(ids) = 1;
  
  % is dq1 small and dq2 relatively large 
  ids = find(abs(dq1)<tol & abs(dq2)>tol );
  thetaL(ids) = 100; % choose theta large
  
  % hopefully the remainder
  ids = find(abs(dq1)>= tol);
  thetaL(ids) = dq2(ids)./dq1(ids);
  
  if(strcmp(method,'minmod'))
      phiL = minmod(1, thetaL);
  elseif(strcmp(method,'MC'))
      phiL = max(0, min(min(0.5*(1+thetaL),2), 2*thetaL));
  elseif(strcmp(method,'vanleer'))
      phiL = (thetaL+abs(thetaL))./(1+abs(thetaL));
  elseif(strcmp(method,'superbee'))
      phiL = max(0, max(min(1,2*thetaL), min(2,thetaL)));
  elseif(strcmp(method,'Fromm'))
      phiL = 0.5*(1+thetaL);
  elseif(strcmp(method,'BeamWarming'))
      phiL = thetaL; % upwind
  elseif(strcmp(method,'LaxWendroff'))
      phiL = ones(size(q)); % downwind
  elseif(strcmp(method, 'upwind'))
      phiL = zeros(size(q));
  elseif(strcmp(method, 'custom'))
     
      phiL = customphi(thetaL);
  else
      disp('WARNING: wrong choice of flux limiter');
  end
  phiR = [phiL(2:N), phiL(1)];
  
  % update cell averages
  q = q - lam*(q-qm1) - 0.5*lam*(1-lam)*(phiR.*(qp1-q) -phiL.*(q-qm1)); 

  if(mod(n,40)==0)
    plot(x,q, '.'); 
    hold on; plot(x-0.5*dx, phiL, 'g-'); hold off;pause(0.02);
  end
end

plot(x, q, 'b.'); 
% have gone one exact period
hold on; plot(x, fluxlimiterexact(x), 'r-'); hold off;
hold on; plot(x-0.5*dx, phiL, 'g-'); hold off;
xlabel('x'); 
legend(sprintf('%s N=%d T=%f', method, N, FinalTime), ...
           'Exact solution', sprintf('flux limiter function: phi_{%s}', method));
axis([-1 1 -1 3])