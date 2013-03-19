% Burgers TVD implementation
% coded by Wenjun @ Nanjing Normal University
% modified by Manuel Diaz @ NTU, 2013.03.19
clc; clear all; close all;

% Parameters
cfl = 1/3;  % Courant Number
tEnd = 2;   % Final time

%% Define Functions
% Flux
f = inline('u.^2/2');
% Superbee Limiter
superbee = inline('max(0,max(min(1,2*r),min(2,r)))');
% Minmod Limiter
minmod = inline('(sign(1)+sign(r))/2.*min(abs(1),abs(r))');

%% Build x-domain
N = 200; x = linspace(0,2,N)'; dx = max(x(2:end)-x(1:end-1));

%% Initial Condition (IC)
u0 = 2+0.5*sin(pi*x); 

%% Main loop
% Initilize Arrays
du1 = zeros(N,1); du2 = zeros(N,1); 
Du1 = zeros(N,1); Du2 = zeros(N,1);
  r = zeros(N,1);  u1 = zeros(N,1);

% Initial time
t = 0; 

% Load IC
u = u0;

while t<tEnd
    % Update iteration time 't',
    dt = cfl*dx/max(abs(u));
    %if t+dt > tEnd; dt = tEnd-t; end;
    t = t+dt; lambda = dt/dx;
    
    %Prepare u^{n}_{i} arrays in space,
    um1 = circshift(u,1); 
    um2 = circshift(um1,1); 
    up1 = circshift(u,-1);
    ubar = 1/2*(u+um1); 
    Lam = lambda*ubar;
    
    % Smoothness function r(x),
    for i = 1:N
        if u(i) == um1(i)
            r(i) = 1;
        elseif ubar(i)>0
            du1(i) = um1(i)-um2(i); du2(i) = u(i)-um1(i);
            r(i) = du1(i)/du2(i);
        elseif ubar<0
            Du1(i) = up1(i)-u(i); Du2(i) = u(i)-um1(i);
            r(i) = Du1(i)/Du2(i);
        end
    end
    
    % Use Limiter:
    %phim = minmod(r); phip = circshift(phim,-1);
    phim = superbee(r); phip = circshift(phim,-1);
    
    Id1 = ubar > 0; Id2 = ubar <= 0;
    if ~isempty(Id1) % u_bar > 0 : Upwind
        u1(Id1) = u(Id1)-lambda*(f(u(Id1))-f(um1(Id1)))-...
            1/2*Lam(Id1).*(1-Lam(Id1)).*(phip(Id1).*(up1(Id1)-u(Id1))-...
            phim(Id1).*(u(Id1)-um1(Id1)));
    end
    if ~isempty(Id2) % u_bar < 0 : Downwind
        u1(Id2) = u(Id2)-lambda*(f(up1(Id2))-f(u(Id2)))-...
            1/2*Lam(Id2).*(1+Lam(Id2)).*(phip(Id2).*(up1(Id2)-u(Id2))-...
            phim(Id2).*(u(Id2)-um1(Id2)));
    end
    
    % Create/Update figure
    plot(x,u1,'o'); axis([x(1) x(end) min(u0)-0.1 max(u0)+0.1]); 
    drawnow;
    
    % Update information
    u = u1;
end