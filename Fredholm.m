%% Solution for the Fredholm integral Equation
% Numerical solution for the equation:
%
% $$f(x)=\phi (x)+\int_{a}^{b}K(x,t)\phi(t)dt$$
%
clear
clc
close all

%% Define range and step zise for dummy variables
% For our case we set x from 0 to pi and t from 0 to 2*pi:
%
hx=0.1;
ht=0.1;
x=0:hx:1*pi;
t=0:ht:1*pi;

%% Define vector f(x_i)
% 
n=length(x);
f=zeros(n,1);
for i=1:n
    f(i)=pi*x(i)^2;
end

%% Define Matriz k(x_i,t_j)

for i=1:n
     for j=1:n
        k(i,j)=0.5*sin(3*x(i))-t(j)*x(i)^2;
     end
end


%% Rectangle Cuadrature
I=eye(n);
G=(I-3*hx*k);
phi=inv(G)*f;
 
%% Exact Solution
%
% $$\phi(x)=sin(3x)$$
%
for i=1:n
    phi_e(i,1)=sin(3*x(i));
end

%% Plot comparison
plot(x,phi,'o',x,phi_e)
