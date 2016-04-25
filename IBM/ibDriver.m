% Driver for 2D Periodic Immeresed Boundary Method Solver
% Solver uses Fourier Space Method of Peskin 
% (see Peskin, Acta Numerica review for method details)
%
%
%  License: This code is free to use for any purposes, provided
%           any publications resulting from the use of this code
%           reference the original code/author.
%
%  Author:  Samuel Isaacson (isaacson@math.utah.edu)
%  Date:    11/2007
%
%  Please notify the author of any bugs, and contribute any
%  modifications or bug fixes back to the original author.
%
%  Disclaimer:
%   This code is provided as is. The author takes no responsibility 
%   for its results or effects.

%
% fluid params:
N   = 64;                % number of mesh points
rho = 1;                 % fluid density
mu  = .01;               % viscosity 
h   = 1 / N;             % mesh width
dt  = .2 * h;            % time step size
Nt  = 2000;               % number of time steps to take

% Lagrangian Params:
dtheta = h / 2;          % Lagrangian mesh spacing 
K  = 1;                  % Stiffness
Kp = .01 / (dt*dt);      % Target point Stiffness

% Euclidean Mesh
X   = (0:(N-1))' * h;    % Coordinate Axis
[X,Y] = meshgrid(X,X);   % X,Y are coordinate matrices 
X = X';
Y = Y';
x = reshape(X, N*N, 1);
y = reshape(Y, N*N, 1);


% initial fluid flow, force:
u0  = zeros(N*N,2); 
f   = zeros(N*N,2);

% initial Lagrangian points:
L    = round( 1 / dtheta );   % number of Lagrangian points
lIdx = (0:(L-1))';
Xl0  = [(.4*cos(2*pi*lIdx*dtheta)+.5) (.15*sin(2*pi*lIdx*dtheta)+.5)]; % point positions


%----------END OF INITIALIZATION CODE-------------%

% declare fluid variables
uX = zeros(N*N,Nt+1);
uY = uX;
u  = u0;
uX(:,1) = u0(:,1);
uY(:,1) = u0(:,2);
p  = zeros(N*N,Nt+1);

% Lagrangian points:
Xlx      = zeros(L, Nt+1);
Xly      = zeros(L, Nt+1);
Xlx(:,1) = Xl0(:,1);
Xly(:,1) = Xl0(:,2);
Xl       = Xl0;

plot(Xl(:,1), Xl(:,2));
hold on;
quiver(reshape(X,N,N), reshape(Y,N,N), reshape(u(:,1),N,N), reshape(u(:,2),N,N) );
hold off;
axis([0 1 0 1]);
pause;

% time step
for(i = 1:Nt)
  fprintf(1, ['i = ' num2str(i) '\n'] );
  
  Xl2    = advanceBoundary(Xl, u, x, y, dt/2, N, h);
  
  F      = calcLagrangianForce(Xl2, L, dtheta, Xl0, K, Kp);
  
  f      = spreadForce(F, Xl2, x, y, dtheta, N, h);    
  
  [u pT] = fluidSolver(u, f, rho, mu, dt, h, N);

  Xl     = advanceBoundary(Xl, u, x, y, dt, N, h);
  
  uX(:,i+1)  = u(:,1);
  uY(:,i+1)  = u(:,2);
  p(:,i+1)   = pT;
  
  Xlx(:,i+1) = Xl(:,1);
  Xly(:,i+1) = Xl(:,2);
  
  
  % plot fluid velocity field and IB points
  plot(Xl(:,1), Xl(:,2));
  hold on;
  quiver(reshape(X,N,N), reshape(Y,N,N), reshape(u(:,1),N,N), reshape(u(:,2),N,N) );
  hold off;
  axis([0 1 0 1]);
  title(['t = ' num2str(dt * i) ', max = ' num2str(max(max(abs(u)))) ', min = ' num2str(min(min(abs(u))))]);
  drawnow
  
end


