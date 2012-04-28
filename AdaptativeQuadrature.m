%% Adaptative Quadrature
% This is an adaptative quadrature algorithm for Simpson and Trapezoidal
% rule.

clear
clc
close all

%% Quadrature grid and initial values
% Integration range and values
%
a=0; % Lower limit
b=1; % Upper limit
n=10; % Initial number of grid divisions
h=(b-a)/n; % Initial step size
x=zeros(1,n+1); % Quadrature points vector
epsilon=1/(b-a)*(1e-8); % Decided error

%% Stablish quadrature point x values
x(1)=a;
for i=1:n
    x(i+1)=x(i)+h;
end

%% Run quadrature function
for i=1:n
    S(i)=SimpsonAQ(x(i),x(i+1),epsilon); %Preseted function of AQ
end
total=sum(S)