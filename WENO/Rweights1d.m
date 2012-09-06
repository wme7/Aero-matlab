function [w] = Rweights1d(u,i,Degree)
% Compute 1D WENO Weights:

% Parameters
%     nx = length(u); 
%      x = 4:nx-3;    
epsilon = 1e-6; % to avoid any alpha's denominator to become zero.

% Initalize Arrays
 beta = zeros(1,Degree);   
alpha = zeros(1,Degree); 
    w = zeros(1,Degree);

% Compute the Weights:
    switch Degree
        
        case{1} % When k = 1; for 2nd order WENO
            
            w = 1; % d = 1
            
        case{2} % When k = 2; for 3rd order WENO
            
            % Smooth Indicators:
            beta(1) = (u(i+1)-u(i))^2; % beta0
            beta(2) = (u(i)-u(i-1))^2; % beta1
            % dr contanst
            d = [2/3 1/3];
            % Smooth coeficients:
            for k = 1:Degree
                alpha(k) = d(k)/(epsilon + beta(k));
            end
            % Weights
            for k = 1:Degree
                w(k) = alpha(k)/(sum(alpha));
            end
            
        case{3} % When k = 3; for 5th order WENO
            
            % Smooth Indicators:
            beta(1) = 13/12*(u(i)-2*u(i+1)+u(i+2))^2 + ...
                1/4*(3*u(i)-4*u(i+1)+u(i+2))^2; % beta(0)
            beta(2) = 13/12*(u(i-1)-2*u(i)+u(i+1))^2 + ...
                1/4*(u(i-1)-u(i+1))^2;          % beta(1)
            beta(3) = 13/12*(u(i-2)-2*u(i-1)+u(i))^2 + ...
                1/4*(3*u(i-2)-4*u(i-1)+u(i))^2; % beta(2)
            % dr contanst
            d = [3/10 6/10 1/10];
            % Smooth coeficients:
            for k = 1:Degree
                alpha(k) = d(k)/(epsilon + beta(k));
            end
            % Weights
            for k = 1:Degree
                w(k) = alpha(k)/(sum(alpha));
            end
            
        otherwise
            error('only options: d = 1, d = 2 and d = 3');
    end