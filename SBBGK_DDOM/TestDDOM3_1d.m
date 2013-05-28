%% DDOM test
% This subroutine was built to test a dynamic integration method
clear all; clc; close all;

%% Choose Statistic:
theta =  0; fprintf('MB statistic \n');
%theta =  1; fprintf('FD statistic \n');
%theta = -1; fprintf('BE statistic \n');

% IC's
%z = 0.3000; u = 2.50; t = 0.1; % random
z = 0.0394; u = 0.50; t = 3.2; % Sod's right
%z = 0.2821; u = 0.75; t = 4.0; % Sod's left

%% Create functions
feq = @(x) 1./((1./z).*exp((x-u).^2./t) + theta); 

% Define particles velocity domain
    % Using GH abscissas as particle velocity domain
    [c,wc] = GaussHermite(10); kc = 1;
        
    % Define Normalized reference
    c_star = (c-u)./sqrt(t);
    
    % GH quadrature function weights:
    wc = wc.*exp((c_star).^2);

% plot comparison of methods
plot(c_star,feq(c),'-ok'); title('DOM over c*')

% 1D  CASE
J = sqrt(t);

%% Apply GH quadrature rule:
fprintf('***************************************************************\n')
fprintf('int(f_eq) \n\n')
tic; n2 = kc*sum(wc'.*feq(c)'.*J); toc; 
fprintf(' int(f_eq*J): %1.12f \n\n',n2); 

fprintf('***************************************************************\n')
fprintf('int(c_star * f_eq) \n\n')
tic; nu2 = kc*sum(wc'.*(c_star+u)'.*feq(c)'.*J^2); toc; 
fprintf(' int(c_star*f_eq*J): %1.12f \n\n',nu2);

fprintf('***************************************************************\n')
fprintf('int(0.5*(c-u)^2 * f_eq) \n\n')
tic; ne2 = 0.5*kc*sum(wc'.*((J*c_star+u).^2)'.*feq(c)'.*J^3); toc; 
fprintf(' int(0.5*(c-u)^2*f_eq*J): %1.12f \n\n',ne2);

fprintf('***************************************************************\n')
fprintf('int(0.5*c^2 * f_eq) \n\n')
tic; nE2 = 0.5*kc*sum(wc'.*((c-u/J).^2)'.*feq(c)'.*J^3); toc; 
fprintf(' int(0.5*c^2*f_eq*J): %1.12f \n\n',nE2);
 
clear c wc n2 feq

%% reference solution
fprintf('***************************************************************\n')
fprintf('Using "quad" function with tolerance = 1E-12\n\n')
tol = 1E-12;
M1 = @(x) 1./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M1,-16,16,tol); toc; 
fprintf('int(f_eq): %1.12f \n\n',n_ref);
fprintf('***************************************************************\n')

M2 = @(x) x./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M2,-16,16,tol); toc; 
fprintf('int(c*f_eq): %1.12f \n\n',n_ref);
fprintf('***************************************************************\n')

M3 = @(x) (0.5*(x-u).^2)./((1/z)*exp((x-u).^2/t) + theta);
tic; n_ref = quad(M3,-16,16,tol); toc; 
fprintf('int(1/2*(c-u)^2*f_eq): %1.12f \n\n',n_ref);

fprintf('***************************************************************\n')

M4 = @(x) (0.5*x.^2)./((1/z)*exp((x-u).^2/t) + theta);
tic; n_ref = quad(M4,-16,16,tol); toc; 
fprintf('int(1/2*c^2*f_eq): %1.12f \n\n',n_ref);

%% Wow!
