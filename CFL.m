%% Homework 10
% We Wish to solve the following partial differential equation:
%
% $\frac{\partial T}{\partial t}+u \frac{\partial T}{\partial x}=\alpha \frac{\partial^2 T}{\partial x^2}$
%
%% Part 1: Pure convection ($\alpha=0$)
%
% Following the condition: $cfl < 1$
%  where $cfl = u \frac{\Delta t}{\Delta x}$

clc
clear all
close all

% Initial Parameters
cfl = 1.1;
u   = 0.08;
dx  = 1/50; % at least 51 grid points
dt  = cfl*dx/u;
t_end  = 8 ;

% grid
x = 0:dx:1;
t = 0:dt:t_end;
X = length(x);
N = length(t);
T = zeros(N,X);

% Initial Condition
T(1,1:12) = 1.-(10.*x(1:12)-1).^2;
T(1,12:X) = 0;

% Plot options
xmin=1; xmax=51; ymin=-0.2; ymax=1.7;

%% Exact Solution
% Computing the Exact solution for comparison purposes:

T_exact0 = zeros(1,X); t_end=0;
for j = 1:X 
    if x(j)-u*t_end <= 0.2 && x(j)-u*t_end >= 0
        T_exact0(j) = 1-(10.*(x(j)-u*t_end)-1)^2;
    else 
        T_exact0(j) = 0;
    end
end

T_exact4 = zeros(1,X); t_end=4;
for j = 1:X 
    if x(j)-u*t_end <= 0.2 && x(j)-u*t_end >= 0
        T_exact4(j) = 1-(10.*(x(j)-u*t_end)-1)^2;
    else 
        T_exact4(j) = 0;
    end
end

T_exact8 = zeros(1,X); t_end=8;
for j = 1:X 
    if x(j)-u*t_end <= 0.2 && x(j)-u*t_end >= 0
        T_exact8(j) = 1-(10.*(x(j)-u*t_end)-1)^2;
    else 
        T_exact8(j) = 0;
    end
end

%% Backward Euler time advancement
% (i) Using Explicit Euler time advancement and second-order central
% difference for the spatial derivate.

% Numerical Scheme:
T(:,1) = 0; T(:,X) = 0;
for n = 1:N-1
    for j = 2:X-1
    T(n+1,j) = T(n,j) - 1/2*cfl*(T(n,j+1) - T(n,j-1));
    end
end

% Plot Figures:
figure
Axis([xmin xmax ymin ymax])
Title(['Backward Euler Scheme',', CFL = ',num2str(cfl)])
hold on
plot(T(1,:),'ro')
plot(T(floor(N/2),:),'b+')
plot(T(N,:),'gs')
plot(T_exact0,'k')
plot(T_exact4,'k')
plot(T_exact8,'k')
legend('@ t = 0','@ t = 4','@ t = 8','Exact Solution')
hold off

%% Leapfrog time advancement
% (ii) Using Leapfrog time advancement and second-order central
% difference for the spatial derivate.
%
% This is a multistep method that is we need to use the first two lines on
% the forward euler to start our new scheeme.

T2 = zeros(N,X);
T2(1,:) = T(1,:);
% Using one sided upwind as our first step
T2(1,1) = 0; T2(1,X) = 0;
for n = 1
    for j = 2:X-1
    T2(n+1,j) = T2(n,j) - cfl*(T2(n,j) - T2(n,j-1));
    end
end

% Numerical Scheme:
T2(:,1) = 0; T2(:,X) = 0;
for n = 2:N-1
    for j = 2:X-1
    T2(n+1,j) = T2(n-1,j) - cfl*(T2(n,j+1) - T2(n,j-1));
    end
end

% Plot Figures:
figure
Axis([xmin xmax ymin ymax])
Title(['Leapfrog Scheme',', CFL = ',num2str(cfl)])
hold on
plot(T2(1,:),'ro')
plot(T2(round(N/2),:),'b+')
plot(T2(N,:),'gs')
plot(T_exact0,'k')
plot(T_exact4,'k')
plot(T_exact8,'k')
legend('@ t = 0','@ t = 4','@ t = 8','Exact Solution')
hold off

%% Lax-Wendroff
% (iii) Using Lax-Wendroff scheme
%
% This is a multistep method that is we need to use the first two lines on
% the forward euler to start our new scheeme.

T3 = zeros(N,X);
T3(1,:) = T(1,:);

% Numerical Scheme:
T3(:,1) = 0; T3(:,X) = 0;
for n = 1:N-1
    for j = 2:X-1
    T3(n+1,j) = T3(n,j) - cfl/2*(T3(n,j+1) - T3(n,j-1))...
        + cfl^2/2*(T3(n,j+1) - 2*T3(n,j) + T3(n,j-1));
    end
end

% Plot Figures:
figure
Axis([xmin xmax ymin ymax])
Title(['Lax-Wendroff Scheme',', CFL = ',num2str(cfl)])
hold on
plot(T3(1,:),'ro')
plot(T3(round(N/2),:),'b+')
plot(T3(N,:),'gs')
plot(T_exact0,'k')
plot(T_exact4,'k')
plot(T_exact8,'k')
legend('@ t = 0','@ t = 4','@ t = 8','Exact Solution')
hold off

%% Assume $u$ was a function of $x$, that is: $u(x)=0.2 sin(\pi x)$ 
% solution:
%
% We know that the maxima and minima of sin(\pi x ) is 1 and -1
% respectively. This implies that the direction of $u$ is changing
% from positive to negative and therefore Backward Euler cannot be used to
% compute with negative velocity $u$. 
% For the leapfrog scheme, however this seems not to be a problem.
% Following the CFL condition we compute the time step as:
%
% $\Delta t < \left | \frac{\Delta x}{0.2 u} \right |$

%% Discussion
% In terms of the stability and accuracy of these schemes 


