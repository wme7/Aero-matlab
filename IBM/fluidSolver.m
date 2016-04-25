function [u, p] = fluidSolver(u0, f, rho, mu, dt, h, N);
% 
% [u] = fluidSolver(u0, f, rho, mu, dt, h, N);
%    
%  Takes one times step of the fluid solver from Dr. Peskin's notes
%
%  Returns:
%     u   = solution at time dt after u0 solution
%     p   = pressure at time dt after u0 solution
%
%  Input:
%     u0  = current solution
%     f   = force at current time + dt/2
%     rho = density
%     mu  = viscosity
%     dt  = time step
%     h   = spatial length
%     N   = number of spatial mesh points in each direction (NxN or NxN
%           grid)
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


persistent D0x D0y L D0xT D0yT LT D0T zeroIdxs isInit;

if( size(isInit) == 0 )
  isInit = 1;
  
  % OPERATORS:
  [D0x, D0y]   = D02DPeriodic(N,h);
  L            = lap2DPeriodic(N,h);
  [D0xT, D0yT] = D02DPeriodicFT(N,h);
  LT           = lap2DPeriodicFT(N,h);
  D0T          = [D0xT D0yT];
  
  % find zero entries in D0T when dotted, replace
  zeroIdxs     = find(dot(D0T,D0T,2) < 2*eps);
  D0T(zeroIdxs,:) = 1;

end


M = N * N;
u = u0;

% 1 - mu*dt/(2*rho) * L
ILT   = 1 - mu*dt / (2*rho) * LT;


% Calulate S(u)u
tmp  = zeros(size(u));
tmp(:,1) = (u(:,1) .* (D0x * u(:,1)) + u(:,2) .* (D0y * u(:,1)) + D0x * (u(:,1) .* u(:,1)) + D0y * (u(:,2) .* u(:,1)) ) / 2;
tmp(:,2) = (u(:,1) .* (D0x * u(:,2)) + u(:,2) .* (D0y * u(:,2)) + D0x * (u(:,1) .* u(:,2)) + D0y * (u(:,2) .* u(:,2)) ) / 2;


% Calculate Transformed rhs
tmp  = u - (dt/2)*tmp + (dt / (2*rho)) * f;
rhs  = fft2D(tmp, N, M);


% Solve for Fourier Space Solutions at t + dt/2:
pT             = 2 * rho / dt * dot( D0T, rhs, 2 ) ./ dot( D0T, D0T, 2 );
pT( zeroIdxs ) = 0;
uT             = ( rhs - (dt/(2*rho)) * [ D0xT .* pT  D0yT .* pT ] ) ./ [ILT ILT]; 


% Transform Back to Real Space
u = ifft2D(uT, N, M);
u = real(u);

% Calculate S(u)u
tmp(:,1) = (u(:,1) .* (D0x * u(:,1)) + u(:,2) .* (D0y * u(:,1)) + D0x * (u(:,1) .* u(:,1)) + D0y * (u(:,2) .* u(:,1)) ) / 2;
tmp(:,2) = (u(:,1) .* (D0x * u(:,2)) + u(:,2) .* (D0y * u(:,2)) + D0x * (u(:,1) .* u(:,2)) + D0y * (u(:,2) .* u(:,2)) ) / 2;


% Calculate Transformed rhs
tmp = u0 - dt*tmp + mu*dt / (2*rho) * (L * u0) + dt/rho * f;
rhs = fft2D(tmp, N, M);


% Solve for Fourier Space Solutions at t + dt:
pT             = rho / dt * dot( D0T, rhs, 2 ) ./ dot( D0T, D0T, 2 );
pT( zeroIdxs ) = 0;
uT             = ( rhs - (dt/(rho)) * [ D0xT .* pT  D0yT .* pT ] ) ./ [ILT ILT]; 


% Transform Back to Real Space
u = ifft2D(uT, N, M);
u = real(u);

p = ifft2D1(pT, N, M);
p = real(p);



function [uT] = fft2D(u, N, M);

uT = [ reshape(fft2( reshape(u(:,1),N,N) ), M, 1) reshape(fft2( reshape(u(:,2),N,N) ), M, 1) ];


function [u] = ifft2D(uT, N, M);

u = [ reshape(ifft2( reshape(uT(:,1),N,N) ), M, 1) reshape(ifft2( reshape(uT(:,2),N,N) ), M, 1) ];


function [u] = ifft2D1(uT, N, M);

u = reshape(ifft2( reshape(uT,N,N) ), M, 1);
