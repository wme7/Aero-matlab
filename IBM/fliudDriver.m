% A driver for just the fluid solver. See ibDriver.m for details...
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


% params:
N   = 64;
rho = 1;
mu  =  .01; %.0005;
h   = 1 / N;
dt  = .75 * h;
Nt  = 500;

% Euclidean Mesh
X   = (0:(N-1))' * h;
[X,Y] = meshgrid(X,X);
X = X';
Y = Y';
x = reshape(X, N*N, 1);
y = reshape(Y, N*N, 1);

% initial fluid flow and force:
u0  = zeros(N*N,2); 
u0 = [(sin(2*pi*x) .* sin(2*pi*y)) (sin(2*pi*x) .* cos(2*pi*y)) ];
f   = zeros(N*N,2);
% u0(30:N:(N*N-N+30),2) = .1;
% u0(31:N:(N*N-N+31),2) = .1;
% u0(32:N:(N*N-N+32),2) = .1;
% u0(33:N:(N*N-N+33),2) = .1;
% u0(257:512,1) = .1;


%----------END OF INITIALIZATION CODE-------------%

% declare fluid variables
uX = zeros(N*N,Nt+1);
uY = uX;
u  = u0;
uX(:,1) = u0(:,1);
uY(:,1) = u0(:,2);
p  = zeros(N*N,Nt+1);

% time step
for(i = 1:Nt)
  i
  
  [u pT] = fluidSolver(u, f, rho, mu, dt, h, N);

  uX(:,i+1) = u(:,1);
  uY(:,i+1) = u(:,2);
  p(:,i+1)  = pT;
end



% plotting stuff:

xmax = max(max(uX));
xmin = min(min(uX));
ymax = max(max(uY));
ymin = min(min(uY));

mmax = max([xmax ymax]);
mmin = min([xmin ymin]);

pmax = max(max(p));
pmin = min(min(p));

hh = figure('DoubleBuffer', 'on');
subplot(2,2,1);
pcolor(X - h/2, Y - h/2, reshape(uX(:,1),N,N));
caxis([mmin mmax]);
colorbar;
title(['u at t = 0']);
xlabel('x');
ylabel('y');

subplot(2,2,2);
pcolor(X - h/2, Y - h/2, reshape(uY(:,1),N,N));
caxis([mmin mmax]);
colorbar;
title(['v at t = 0']);
xlabel('x');
ylabel('y');


subplot(2,2,3);
pcolor(X - h/2, Y - h/2, reshape(p(:,1),N,N));
caxis([pmin pmax]);
colorbar;
title(['p at t = 0']);
xlabel('x');
ylabel('y');

drawnow 
pause

for(i = 2:(Nt+1))
  
  subplot(2,2,1);
  pcolor(X - h/2, Y - h/2, reshape(uX(:,i),N,N));
  caxis([mmin mmax]);
  colorbar;
  title(['u at t = ' num2str(dt*(i-1))]);
  xlabel('x');
  ylabel('y');
  
  subplot(2,2,2);
  pcolor(X - h/2, Y - h/2, reshape(uY(:,i),N,N));
  caxis([mmin mmax]);
  colorbar;
  title(['v at t = ' num2str(dt*(i-1))]);
  xlabel('x');
  ylabel('y');

  subplot(2,2,3);
  pcolor(X - h/2, Y - h/2, reshape(p(:,i),N,N));
  caxis([pmin pmax]);
  colorbar;
  title(['p at t = ' num2str(dt*(i-1))]);
  xlabel('x');
  ylabel('y');
  
  
  drawnow;
  pause(.1);
end