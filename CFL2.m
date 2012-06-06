%% Homework 10
% We Wish to solve the following partial differential equation:
%
% $\frac{\partial T}{\partial t}+u \frac{\partial T}{\partial x}=\alpha \frac{\partial^2 T}{\partial x^2}$
%
%% Part 2: Convection and Diffusion (Assume $\alpha=0.001$)
% Now using the same initial boundary conditions in Part 1, Solve the
% convection-diffusion equation. Repeat part(a): i and ii with the
% addition of the second-order central difference for the diffusion term.
%
% Following the condition: $cfl < 1$
%  where $cfl = u \frac{\Delta t}{\Delta x}$
%
% And now define $\alpha*\frac{\Delta t}{(\Delta x)^2}$ in the range [0,1]

clc
clear all
close all

% Initial Parameters
cfl    = 0.7;
u      = 0.08;
dx     = 1/50; % at least 51 grid points
dt     = cfl*dx/u;
alpha  = 0.001;
beta   = alpha*dt/dx^2;
t_end  = 8;

% Grid
x = 0:dx:1;
t = 0:dt:t_end;
X = length(x);
N = length(t);
T = zeros(N,X);

% Initial Condition
T(1,1:12)  = 1.-(10.*x(1:12)-1).^2;
T(1,12:51) = 0;

% Plot options
xmin=1; xmax=51; ymin=-0.2; ymax=1.7;

%% Barckward Euler time advancement
% a) Using Explicit Euler time advancement and second-order central
% difference for the spatial derivate.

% Numerical Scheme:
T(:,1) = 0; T(:,X) = 0;
for n = 1:N-1
    for j = 2:X-1
    T(n+1,j) = T(n,j) - 1/2*cfl*(T(n,j+1) - T(n,j-1))...
        + beta*(T(n,j+1) - 2*T(n,j) + T(n,j-1));
    end
end

% Plot Figures:
figure
Axis([xmin xmax ymin ymax])
Title(['Backward Euler Scheme',', CFL = ',num2str(cfl),', Beta = ',num2str(beta)])
hold on
plot(T(1,:),'ro')
plot(T(floor(N/2),:),'b+')
plot(T(N,:),'gs')
legend('@ t = 0','@ t = 4','@ t = 8')
hold off

%% Leapfrog time advancement
% b) Using Leapfrog time advancement and second-order central
% difference for the spatial derivate and a central difference for the
% diffusion term.
%
% This is a multistep method that is we need to use the first two lines on
% the forward euler to start our new scheme.

T2 = zeros(N,X);
T2(1,:) = T(1,:);
T2(2,:) = T(2,:);
T2(1,1) = 0; T2(1,X) = 0;

% Numerical Scheme:
T2(:,1) = 0; T2(:,X) = 0;
for n = 2:N-1
    for j = 2:X-1
    T2(n+1,j) = T2(n-1,j) - cfl*(T2(n,j+1) - T2(n,j-1))...
        + 2*beta*(T2(n,j+1) - 2*T2(n,j) + T2(n,j-1));
    end
end

% Plot Figures:
figure
%Axis([xmin xmax ymin ymax])
Title(['Leapfrog Scheme',', CFL = ',num2str(cfl),', Beta = ',num2str(beta)])
hold on
plot(T2(1,:),'-.ro')
plot(T2(floor(N/2),:),'-.b+')
plot(T2(N,:),'-.gs')
legend('@ t = 0','@ t = 4','@ t = 8')
hold off

%% Leapfrog time advancement with the diffusion term lagged in time.
% c) Using Leapfrog time advancement and second-order central
% difference for the spatial derivate and a central difference for the
% diffusion term. But this time the diffusion term would be evaluated at
% step "n-1" rather than "n".
%
% This is a multistep method that is we need to use the first two lines on
% the forward euler to start our new scheme.

T3 = zeros(N,X);
T3(1,:) = T(1,:);
T3(2,:) = T(2,:);
T3(1,1) = 0; T3(1,X) = 0;

% Numerical Scheme:
T3(:,1) = 0; T3(:,X) = 0;
for n = 2:N-1
    for j = 2:X-1
    T3(n+1,j) = T3(n-1,j) - cfl*(T3(n,j+1) - T3(n,j-1))...
        + 2*beta*(T3(n-1,j+1) - 2*T3(n-1,j) + T3(n-1,j-1));
    end
end

% Plot Figures:
figure
Axis([xmin xmax ymin ymax])
Title(['Leapfrog with diff-term lagged Scheme',', CFL = ',num2str(cfl),', Beta = ',num2str(beta)])
hold on
plot(T3(1,:),'ro')
plot(T3(floor(N/2),:),'b+')
plot(T3(N,:),'gs')
legend('@ t = 0','@ t = 4','@ t = 8')
hold off

%% Discusion
% How the presence of the diffusion term affected the physical behavior of
% the solution and stability properties of the numerical solution?

