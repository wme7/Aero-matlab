%% Explicit Euler
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
% initial condition (v(0))
v = zeros(n,1);
alpha = zeros(n,1);
beta  = zeros(n,1);

v(1)=1; 
for i=1:n-1
    % compute variables alpha(t) and beta(t)
    alpha(i) = 3*t(i)/(1+t(i));
    beta(i)  = 2*(1+t(i))^3*exp(-t(i));
    % main discrete equation
    v(i+1) = v(i) + h*(-alpha(i)*v(i) + beta(i));
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