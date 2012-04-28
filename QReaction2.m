%% Reaction Solver
% Homework 7, Problem 2
clc
close all
clear all
%% Convencion
% c(1) = Cs; c(2) = Ce, c(3) = Ces, c(4) = Cp

%% Solving with a linearized trapezoidal method:
% Discretizing time:
ti = 0;   % initial time
tf = 1E2; % final time
h = 0.1;  % time step
t = ti:h:tf;

% Discretizing concentration variables:
n = length(t); c = zeros(n,4); % i-concentration


%% Initial Conditions
c(1,:) = [1 5.0e-5 0.0 0.0];
k      = [2.0E3 10E-3 10.0];

%% Crank Nicolson Scheeme
for i=1:n-1
    % f functions
    f1 = -k(1)*c(i,1)*c(i,2) +  k(2)*c(i,3);
    f2 = -k(1)*c(i,1)*c(i,2) + (k(2)-k(3))*c(i,3);
    f3 =  k(1)*c(i,1)*c(i,2) - (k(2)+k(3))*c(i,3);
    f4 =  k(3)*c(i,3);
    f = [f1 f2 f3 f4]';
    
    % Jacobian Matrix
    A = [-k(1)*c(i,2)  -k(1)*c(i,1)  k(2)       0;...
         -k(1)*c(i,2)  -k(1)*c(i,1)  k(2)-k(3)  0;...
          k(1)*c(i,2)  -k(1)*c(i,1) -k(2)+k(3)  0;...
          0            0             k(3)       0];
         
    c(i+1,:) = (c(i,:)'+ h/2*inv((eye(4) - h/2*A))*(f + f))';
    %c(i+1,:) = (c(i,:)'+ h/2*((eye(4) - h/2*A))\(f + f))';
end

%% Make nice figure

loglog(t,c(:,1),'blue',t,c(:,2),'red',t,c(:,3),'green',t,c(:,4),'black');
axis([ti tf 10e-5 10e6]);
grid on
legend('Cs','Ce','Ces','Cp');


%     % Jacobian Matrix
%     dc1 = c(i+1,1)-c(i,1);
%     dc2 = c(i+1,2)-c(i,2);
%     dc3 = c(i+1,3)-c(i,3);
%     dc4 = c(i+1,4)-c(i,4);
%     df1 = f(i+1,1)-f(i,1);
%     df2 = f(i+1,2)-f(i,2);
%     df3 = f(i+1,3)-f(i,3);
%     df4 = f(i+1,4)-f(i,4);
% 
%     A = [df1/dc1 df1/dc2 df1/dc3 df1/dc4;...
%          df2/dc1 df2/dc2 df2/dc3 df2/dc4;...
%          df3/dc1 df3/dc2 df3/dc3 df3/dc4;...
%          df4/dc1 df4/dc2 df4/dc3 df4/dc4];