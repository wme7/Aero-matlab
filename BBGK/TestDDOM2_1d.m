%% DDOM test
% This subroutine was abuilt to test a dynamic integration method
clear all; clc; close all;

%% Choose Statistic:
%theta =  0; fprintf('MB statistic \n');
%theta =  1; fprintf('FD statistic \n');
theta = -1; fprintf('BE statistic \n');

% IC's
%z = 0.3000; u = 2.50; t = 0.1; % random
%z = 0.0394; u = 0.00; t = 3.2; % Sod's right
z = 0.2821; u = -1.35; t = 4.0; % Sod's left

%% Create functions
feq = @(x) 1./((1./z).*exp((x-u).^2./t) + theta); 

% Define particles velocity domain
quadrature = 2;
switch quadrature
    case{1} % Newton Cotes
    % Uniform particle velocity domain
    [v,wv,kv] = cotes_xw(-20,20,100,5);
    [c,wc,kc] = cotes_xw(-20,20,100,5);
    
    % Define Normalized reference
    c_star = (c-u)./sqrt(t);
    
    case{2} % Gauss-Hermite
    % Using GH abscissas as particle velocity domain
    [v,wv] = GaussHermite(20); kv = 1;
    [c,wc] = GaussHermite(5); kc = 1;
        
    % Define Normalized reference
    c_star = (c-u)./sqrt(t);
    
    % GH quadrature function weights:
    wv = wv.*exp(v.^2); 
    wc = wc.*exp((c_star).^2);
    %wc = wc.*exp((c_star).^2 + log(theta));

end

% Define normalized domain
fprintf('max c_star: %1.2f \n',max(c_star))
fprintf('min c_star: %1.2f \n',min(c_star))

% plot comparison of methods
subplot(1,2,1); plot(v,feq(v),'-sr'); title('Traditional DOM over c')
subplot(1,2,2); plot(c_star,feq(c),'-ok'); title('DOM over c*')

% 1D  CASE
J = sqrt(t);

%% Apply GH quadrature rule:
fprintf('***************************************************************\n')
fprintf('int(f_eq) \n\n')
tic; n1 = kv*sum(wv'.*feq(v)'); toc; 
tic; n2 = kc*sum(wc'.*feq(c)'.*J); toc; 
fprintf(' int(f_eq): %1.12f \n',n1);
fprintf(' int(f_eq*J): %1.12f \n\n',n2); 

% fprintf('***************************************************************\n')
% fprintf('int(c*f_eq) \n\n')
% tic; n3 = sum(v'.*wv'.*feq'); toc; 
% tic; n4 = sum(c'.*wc'.*feqc'.*J); toc; 
% fprintf('\n int(c*f_eq): %1.12f \n',n3);
% fprintf(' int(c*f_eqc*J): %1.12f \n\n',n4);
% 
% fprintf('***************************************************************\n')
% fprintf('int(1/2*(c-u)^2*f_eq) \n\n')
% tic; n5 = sum(0.5*(v-u).^2'.*wv'.*feq'); toc; 
% tic; n6 = sum(0.5*(c-u).^2'.*wc'.*feqc'.*J); toc; 
% fprintf('\n int(1/2*(c-u)^2*f_eq): %1.12f \n',n5);
% fprintf(' int(1/2*(c-u)^2*f_eqc*J): %1.12f \n\n',n6);
% 
% fprintf('***************************************************************\n')
% fprintf('int(1/2*c^2*f_eq) \n\n')
% tic; n7 = sum(0.5*v.^2'.*wv'.*feq'); toc; 
% tic; n8 = sum(0.5*c.^2'.*wc'.*feqc'.*J); toc; 
% fprintf('\n int(1/2*c^2*f_eq): %1.12f \n',n7);
% fprintf(' int(1/2*c^2*f_eqc*J): %1.12f \n\n',n8);
% 
%% reference solution
fprintf('***************************************************************\n')
fprintf('Using "quad" function with tolerance = 1E-12\n\n')
tol = 1E-12;
M1 = @(x) 1./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M1,-10,10,tol); toc; 
fprintf('int(f_eq): %1.12f \n\n',n_ref);
% M2 = @(x) x./((1/z)*exp((x-u).^2/t) + theta); 
% tic; n_ref = quad(M2,-10,10,tol); toc; 
% fprintf('int(c*f_eq): %1.12f \n\n',n_ref);
% M3 = @(x) (0.5*(x-u).^2)./((1/z)*exp((x-u).^2/t) + theta);
% tic; n_ref = quad(M3,-10,10,tol); toc; 
% fprintf('int(1/2*(c-u)^2*f_eq): %1.12f \n\n',n_ref);
% M4 = @(x) (0.5*x.^2)./((1/z)*exp((x-u).^2/t) + theta);
% tic; n_ref = quad(M4,-10,10,tol); toc; 
% fprintf('int(1/2*(c-u)^2*f_eq): %1.12f \n\n',n_ref);
fprintf('***************************************************************\n')
% 
% %% Macroscopic conditions recovered
% fprintf('\n\n');
% fprintf('rho = %1.8f',n2);
% fprintf(' ux = %1.8f',n4/n2);
% fprintf('  t = %1.8f',4*n8/n2 - 2*(n4/n2)^2); % t = 4*E./n - 2*u.^2;
% fprintf('  z = %1.8f',n2/sqrt(pi*(4*n8/n2 - 2*(n4/n2)^2))); % r = n./sqrt(pi.*t);
% fprintf('\n\n');

% WOW!