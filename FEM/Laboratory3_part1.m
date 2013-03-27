%% Laboratory No.3, Part 1: Gauss Integration
% Coded by Manuel Diaz, 2013.03.20
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute and plot integration errors with 1 to 5 Guass quadrature points

% Prepare figure 1
figure(1)
for j = 1:2
    % Choose polynomial
    pol = j;
    switch pol
        case{1} % polynomial #1
            p{1} = 'p(x) = 20+x+2*x.^2+3*x.^3+4*x.^4+5*x.^5+6*x.^6+7*x.^7';
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
    name=['Integrating with N-Gauss quadrature points: ',p(j)];title(name);
    ylabel('Error %'); xlabel('Number of gauss points')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
