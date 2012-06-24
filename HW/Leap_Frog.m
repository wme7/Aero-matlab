%% Leapfrog vs. Adams Bashforth
% Routine for folving the pendulum problem.
% 
% for solving 
%
% $\frac{d \theta}{dt}=-\frac{g}{l} \theta$
%
% Assuming for the period 0<t<6 with initial conditions:
%
% $\theta(0) = \pi/18$ and $\theta'(0) = 0$
clc
clear all
close all
%% Constants
g = 9.81; % m/s^2
l = 0.60; % m

%% Discretization of t
ti = 0;   % Initial Time
tf = 6;   % Final Time

% Time step, evaluate for h= 0.1,0.2,0.5
h = 0.5; 

% Discretization of t
t = ti:h:tf;
n = length(t);

%% Computing Theta
theta1 = zeros(n,1);
theta2 = zeros(n,1);
theta1(1) = pi/18; %initial condition
theta2(1) = 0;     %initial condition
% explicit euler as our first step
theta1(2) = theta1(1) + h*theta2(1);
theta2(2) = theta2(1) + h*(-g/l)*theta1(1);

for i=2:n-1
    %Leap-Frog scheeme
    theta1(i+1) = theta1(i-1) + 2*h*theta2(i);
    theta2(i+1) = theta2(i-1) + 2*h*(-g/l)*theta1(i);
end

theta11 = zeros(n,1);
theta22 = zeros(n,1);
theta11(1) = pi/18; %initial condition
theta22(1) = 0;     %initial condition
% explicit euler as our first step
theta11(2) = theta11(1) + h*theta22(1);
theta22(2) = theta22(1) + h*(-g/l)*theta11(1);

for i=2:n-1
    %Adams-Bashforth scheeme
    theta11(i+1) = theta11(i) + 3/2*h*theta22(i)...
        - h/2*theta22(i-1);
    theta22(i+1) = theta22(i) + 3/2*h*(-g/l)*theta11(i)... 
        - h/2*(-g/l)*theta11(i-1);
end
%% Exact solution
theta = zeros(n,1);
theta = theta1(1)*cos( sqrt(g/l).*t );

%% Make figures
hold on
grid on
axis([0 6 -0.4 0.4]);
plot(t,theta,'blue');
plot(t,theta1,'red');
plot(t,theta11,'green');
legend('Exact','Leapfrog','Adams Bashforth');
xlabel('t'); ylabel('Theta(t)')
grid on
title(h)
hold off