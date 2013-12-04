function [dfdci] = GaussSlope(ci,fi)
% This function is used to compute the slope values at ci poinst give the
% values fi of a gaussian-like distribution function.
% Requirement: ci, fi and dfdci must be arryas of the same size.

%% Arrays information
[nv,nx] = size(fi);

% Initialize variable
dfdci = zeros(nv,nx);

%% Fitting
clog = ci; flog = log(fi); 
for i = 1:nx
    % Our final function as computed in classical statistics,
    % y = A*exp(-(x-mu).^2/(2*sigma^2));
    
    p = polyfit(clog(:,i),flog(:,i),2);
    A2 = p(1); A1 = p(2); A0 = p(3);
    sigma = sqrt(-1/(2*A2));
    mu = A1*sigma^2;
    A = exp(A0+mu^2/(2*sigma^2));
    
    %% Derivative
    dfdci(:,i) = A*(mu-ci(:,i))/sigma^2.*exp(-(ci(:,i)-mu).^2/(2*sigma^2));
    
end