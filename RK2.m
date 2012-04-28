%% Runge-Kutta 2
% 
% for solving 
%
% $\frac{dv}{dt}=-\alpha (t)v+\beta (t)$
%
% where $\alpha (t) = \frac{3t}{1-t}$ and $\beta (t)=2(1-t)^3 e^{-t}$
%
% Assuming v(0)=1.0 for the period 0<t<15

clc
clear 
%close all

%% Discretization of t
ti=0;    % Initial Time
tf=150;   % Final Time

% Time step 
% h = 0.2; 
% h = 0.8;
% h = 1.1;
 h = 0.68;

% Discretization of t
i=1;
t(1)=ti;
while t(i)<tf  
        i=i+1;
        t(i)=t(i-1)+h;
end
t(i)=tf; clear i;
n=length(t);

%% Discretization of v

v = zeros(n,1);
alpha = zeros(n,1);
beta  = zeros(n,1);
t2 = zeros(n,1);
alpha2 = zeros(n,1);
beta2  = zeros(n,1);

for i=1:n-1
    t2(i)     = t(i)+1/2*h;
    % compute variables alpha(t) and beta(t)
    alpha(i)  = 3*t(i)/(1+t(i));
    alpha2(i) = 3*t2(i)/(1+t2(i)); % alpha(t + 1/2h)
    beta(i)   = 2*(1+t(i))^3*exp(-t(i));
    beta2(i)  = 2*(1+t2(i))^3*exp(-t2(i)); % beta(t + 1/2h)
end

% initial condition (v(0))
v(1)=1; 

for i=1:n-1
    % main discrete equation
    k1 = h*(-alpha(i)*v(i)+beta(i));
    k2 = h*(-alpha2(i)*(v(i)+1/2*k1)+beta2(i));
    v(i+1) = v(i)+k2;
end

%% Exact solution
v2 = zeros(n,1);
for i=1:n
     v2(i) = exp(-t(i))*(1+t(i))^3;
end
error=v2-v;
figure
plot(t,v,'.',t,v2,'-')
%figure
%plot(t,error,'.')