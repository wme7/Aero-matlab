%% Test Newton Method 2
% Test Newton_n_sytem subroutine

% Clear memory
clear all; close all; clc;

% display format
format long

%% Examples
example = 2;
switch example
    case{1} % 2-by-2 Polynomials,
        % Define vector of functions 'f' and matrix of partial derivates 'df',
        f = inline('[x(1)^3+x(2)-1 ; x(2)^3-x(1)+1]');
        df = inline('[3*x(1)^2, 1  ; -1,  3*x(2)^2]');

        % Plot figures,
        g1 = 'x^3+y-1'; g2 = 'y^3-x+1'; xy_range = [-4,4,-4,4];
        figure(1); grid on;
        hold on;
        ez1 = ezplot(g1,xy_range);
        ez2 = ezplot(g2,xy_range); 
        hold off;
        legend('g1(x,y)=x^3+y-1','g2(x,y)=y^3-x+1'); 
        set(ez1,'color',[0 0 1]);
        set(ez2,'color',[1 0 0]);
        title('Solving x^3+y-1 & y^3-x+1 for {x,y}')
        
        % Initial Guess,
        x0 = [.5;.5] ;
        
    case{2} % 2-by-2 Nonlinear system,
        % Define vector of functions 'f' and matrix of partial derivates 'df',
        f = inline('[10*x(1)^2+sin(x(2))-20 ; x(1)^4+5*x(2)-6]');
        df = inline('[20*x(1), cos(x(2))  ; 4*x(1)^3, 5]');
        
        % Plot figures,
        g1 = '10*x^2+sin(y)-20'; g2 = 'x^4+5*y-6'; xy_range = [-8,8,-8,8];
        figure(1); grid on;
        hold on;
        ez1 = ezplot(g1,xy_range);
        ez2 = ezplot(g2,xy_range); 
        hold off;
        legend('g1(x,y)=10x^2+sin(y)-20','g2(x,y)=x^4+5y-6',1); 
        set(ez1,'color',[0 0 1]);
        set(ez2,'color',[1 0 0]);
        title('Solving 10x^2+sin(y)=20 & x^4+5y=6 for {x,y}')
        
        % Initial Guess,
        x0 = [1;1];
        
        case{3} % 3-by-3 Polynomial System,
        % Define vector of functions 'f' and matrix of partial derivates 'df',
        f = inline('[x(1)^3-2*x(2)-2 ; x(1)^3-5*x(3)^2+7 ; x(2)*x(3)^2-1]');
        df = inline('[3*x(1)^2 -2 0;3*x(1)^2 0 -10*x(3);0 x(3)^2 2*x(2)*x(3)]');
      
        % Initial Guess,
        x0 = [1;1;1];
end

%% Convergence criteria
    % Define max number of iterations,
    n_max = 8;

    % Define tolerance,
    tol = 1e-6;
    tol = tol*ones(size(2,1));

%% Use Newton_n_sytem.m func,
tic;
x = x0; x_next = Newton_n_system(f,df,x); n = 1;
while and((n < n_max),(abs(x_next - x) > tol))
    x = x_next;
    x_next = Newton_n_system(f,df,x);
    n = n + 1;
end
fprintf('My solver: %1.12f \n',x) % show me the result
t1 = toc; % time for Costumb solver

%% Using Matlab solver,
tic;
x = x0; x = fsolve(f,x);
fprintf('Matlab solver: %1.12f \n',x) % show me the result
t2 = toc; % time for Matlab solver

% how fast is our Newton Solver,
fprintf('Costum Solver is: %1.1f times faster than Matlab solver \n',t2/t1)