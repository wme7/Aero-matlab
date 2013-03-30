%% DDOM test
% This subroutine was abuilt to test a dynamic integration method
clear all; clc; close all;

%% Choose Statistic:
theta =  0; fprintf('MB statistic \n');
%theta =  1; fprintf('FD statistic \n');
%theta = -1; fprintf('BE statistic \n');

% IC's
z = 0.3000; u = 2.50; t = 0.1; % random
%z = 0.0394; u = 0.00; t = 3.2; % Sod's right
%z = 0.2821; u = -1.35; t = 4.0; % Sod's left

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

clear c wc n2 feq

%% reference solution
fprintf('***************************************************************\n')
fprintf('Using "quad" function with tolerance = 1E-12\n\n')
tol = 1E-12;
M1 = @(x) 1./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M1,-16,16,tol); toc; 
fprintf('int(f_eq): %1.12f \n\n',n_ref);

%% Wow
