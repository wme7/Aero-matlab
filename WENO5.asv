function v = weno5(f,a,b,N,nodes)
%% Weighted Essentially Non-Oscilatory 5
% WENO 5 implementation subroutine for computing fluxes:
% To solve for a 1D scalar advection equation. Algorith based on:
%
%   Chi-Wang Shu; High-Order ENO and WENO schemes for Computational
%   Fluid Dynamics, High-Order Methods for Computational Physics.
%   Springer 1999.

v = f(x);  k = 2;  epsilon = 1e-25;

i=1;

%% Compute the weight:
switch k 
    
    case{1} % When k = 1; for second order WENO
    d0 = 1;
    w0 = 1; % not sure yet
        
    case{2} % When k = 2; for third order WENO
        % Smooth Indicators:
    beta0= (v(i+1)-v(i))^2; 
    beta1= (v(i)-v(i-1))^2;
        % dr contanst
    d0 = 2/3; d1 = 1/3;
        % Smooth coeficients:
    alpha0 = d0/(epsilon + beta0);
    alpha1 = d1/(epsilon + beta1);
        % Weights 
    w0 = alpha0/(alpha0 + alpha1);
    w1 = alpha1/(alpha0 + alpha1);
    
    case{3} % When k = 3; for fifth order WENO
        % Smooth Indicators:
    beta0= 13/12*(v(i)-2*v(i+1)+vi(i+2))^2 + 1/4*(3*v(i)-4*v(i+1)+v(i+2))^2;
    beta1= 13/12*(v(i-1)-2*v(i)+vi(i+1))^2 + 1/4*(3*v(i)-4*v(i+1)+v(i+2))^2;
    beta2= 13/12*(v(i-1)-2*v(i-1)+vi(i))^2 + 1/4*(3*v(i)-4*v(i+1)+v(i+2))^2;
        % dr contanst
    d0 = 3/10; d1 = 3/5; d2 = 1/10; 
        % Smooth coeficients:
    alpha0 = d0/(epsilon + beta0);
    alpha1 = d1/(epsilon + beta1);
    alpha2 = d2/(epsilon + beta2);
        % Weights 
    w0 = alpha0/(alpha0 + alpha1 + alpha2);
    w1 = alpha1/(alpha0 + alpha1 + alpha2);
    w2 = alpha2/(alpha0 + alpha1 + alpha2);
    
    otherwise
        error('only options: 1, 2 and 3');
end



