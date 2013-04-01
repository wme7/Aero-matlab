%% Laboratory No.3, Part 1: Gauss Integration
% Coded by Manuel Diaz, 2013.03.20
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute and plot integration errors with 1 to 5 Guass quadrature points

% Prepare figure 1
figure(1)
for j = 1:2
    % Choose polynomial
    pol = j;
    switch pol
        case{1} % polynomial #1
            p{1} = 'p(x) = 20+x+2x^2+3x^3+4x^4+5x^5+6x^6+7x^7';
            P = @(x) 20+x+2*x.^2+3*x.^3+4*x.^4+5*x.^5+6*x.^6+7*x.^7;
        case{2} % polynomial #2
            p{2} = '1./(1+x.^2)';
            P = @(x) 1./(1+x.^2);
    end
    
    % define interval of integration
    a = 1;
    b = 2*sqrt(3);
    
    % perform gauss integration
    ngp = 1:5; % = {1,2,3,4,5}
    approx = zeros(1,ngp);
    for i = ngp
        approx(i) = gaussint(P,a,b,i);
    end
    
    % compute exact solution
    syms x;
    exact = subs(int(P,x,a,b));
    
    % compute and plot error
    subplot(1,2,j)
    error = abs((exact-approx)/exact)*100;
    plot(ngp,error,'ob','MarkerFaceColor','g','MarkerSize',10);
    name=['Integrating with N-Gauss quad. points: ',p(j)];title(name);
    ylabel('Error %'); xlabel('Number of gauss points')
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve the loaded Tapered beam
% Solve it's displacement and stress and compare with it exact solution.

% Load cases per number of elements
elements = [1 2 4 8 16];

% Solve and plot Displacement
figure(2)
hold on
for j = 1:5
   [D,S,x] = Driver1d(elements(j));
   markers = {':+r','--og','-.*b',':xm','--sk','-.dw'}; 
   plot(x,D,markers{j})
end
title('Displacement vs x'); ylabel('displacment (m)'); xlabel('x (m)');

% Compute exact solution
exact_D= 1/35000*log(x+0.5);
plot(x,exact_D,'-r')
legend('1 element','2 element','4 element','8 element','16 element','exact',2)
hold off

% Compute and plot stress
figure(3)
hold on
for j = 1:5
   [D,S,x] = Driver1d(elements(j));
   markers = {':+r','--og','-.*b',':xm','--sk','-.dw'}; 
   S = [S;S(end)]; % repeat final value to match number of nodes
   stairs(x,S,markers{j})
   clear S
end
title('Stress vs x'); ylabel('displacment (m)'); xlabel('x (m)');

% Compute exact solution
exact_S = 1/35000*70e9./(x+0.5);
plot(x,exact_S,'-r')
legend('1 element','2 element','4 element','8 element','16 element','exact',1)
hold off
