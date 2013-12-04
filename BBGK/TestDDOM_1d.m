%% DDOM test
% This subroutine was abuilt to test a dynamic integration method
clear all; clc; close all;

% Statistic:
theta =  0; fprintf('MB statistic');
%theta =  1; fprintf('FD statistic');
%theta = -1; fprintf('BE statistic');

% IC's
%z = 0.3000; u = 1.00; t = 0.1; % random
%z = 0.0394; u = 0.00; t = 3.2; % Sod's right
z = 0.2821; u = 0.75; t = 4.0; % Sod's left

% Create functions
M = @(c) 1./((1./z).*exp((c-u).^2./t) + theta); 
Mc = @(c_star) 1./((1./z).*exp(c_star.^2) + theta);

% Weight and C* values
nv = input('Number of quadrature desided quadrature points: ');
%nv = 3; 
fprintf('quadrature points used: %1.0f\n\n',nv);
[c_star,wc] = GaussHermite(nv);
[v,wv] = GaussHermite(nv);

%% 1D  CASE
J = sqrt(t);

% Transformation rule of c and C*
c = sqrt(t) * c_star + u; 

% Evaluate

feq  = M(v) %#ok<NOPTS>
feqc = Mc(c_star) %#ok<NOPTS>

%if feq == feqc; fprintf('feq & feqc agree \n\n'); else fprintf('feq & feqc disgree \n\n'); end;

% GH quadrature function weight:
wv = wv.*exp(v.^2); 
wc = wc.*exp(c_star.^2);
fprintf('\n');

% Apply GH quadrature rule:
fprintf('***************************************************************\n')
fprintf('int(f_eq) \n\n')
tic; n1 = sum(wv'.*feq'); toc; 
tic; n2 = sum(wc'.*feqc'.*J); toc; 

fprintf('\n int(f_eq): %1.12f \n',n1);
fprintf(' int(f_eqc*J): %1.12f \n\n',n2); 

fprintf('***************************************************************\n')
fprintf('int(c*f_eq) \n\n')
tic; n3 = sum(v'.*wv'.*feq'); toc; 
tic; n4 = sum(c'.*wc'.*feqc'.*J); toc; 
fprintf('\n int(c*f_eq): %1.12f \n',n3);
fprintf(' int(c*f_eqc*J): %1.12f \n\n',n4);

fprintf('***************************************************************\n')
fprintf('int(1/2*(c-u)^2*f_eq) \n\n')
tic; n5 = sum(0.5*(v-u).^2'.*wv'.*feq'); toc; 
tic; n6 = sum(0.5*(c-u).^2'.*wc'.*feqc'.*J); toc; 
fprintf('\n int(1/2*(c-u)^2*f_eq): %1.12f \n',n5);
fprintf(' int(1/2*(c-u)^2*f_eqc*J): %1.12f \n\n',n6);

fprintf('***************************************************************\n')
fprintf('int(1/2*c^2*f_eq) \n\n')
tic; n7 = sum(0.5*v.^2'.*wv'.*feq'); toc; 
tic; n8 = sum(0.5*c.^2'.*wc'.*feqc'.*J); toc; 
fprintf('\n int(1/2*c^2*f_eq): %1.12f \n',n7);
fprintf(' int(1/2*c^2*f_eqc*J): %1.12f \n\n',n8);

%% reference solution
fprintf('***************************************************************\n')
fprintf('Using "quad" function with tolerance = 1E-12\n\n')
tol = 1E-12;
M1 = @(x) 1./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M1,-10,10,tol); toc; 
fprintf('int(f_eq): %1.12f \n\n',n_ref);
M2 = @(x) x./((1/z)*exp((x-u).^2/t) + theta); 
tic; n_ref = quad(M2,-10,10,tol); toc; 
fprintf('int(c*f_eq): %1.12f \n\n',n_ref);
M3 = @(x) (0.5*(x-u).^2)./((1/z)*exp((x-u).^2/t) + theta);
tic; n_ref = quad(M3,-10,10,tol); toc; 
fprintf('int(1/2*(c-u)^2*f_eq): %1.12f \n\n',n_ref);
M4 = @(x) (0.5*x.^2)./((1/z)*exp((x-u).^2/t) + theta);
tic; n_ref = quad(M4,-10,10,tol); toc; 
fprintf('int(1/2*(c-u)^2*f_eq): %1.12f \n\n',n_ref);
fprintf('***************************************************************\n')

%% Macroscopic conditions recovered
fprintf('\n\n');
fprintf('rho = %1.8f',n2);
fprintf(' ux = %1.8f',n4/n2);
fprintf('  t = %1.8f',4*n8/n2 - 2*(n4/n2)^2); % t = 4*E./n - 2*u.^2;
fprintf('  z = %1.8f',n2/sqrt(pi*(4*n8/n2 - 2*(n4/n2)^2))); % r = n./sqrt(pi.*t);
fprintf('\n\n');

%% plot comparison solution
hold on
x = -7:0.1:7; grid on; 
plot(x,M(x),'--r',...
            'LineWidth',2); 
plot(x,Mc(x),'-.k',...
            'LineWidth',2);
plot(v,M(v),'ro',...
            'MarkerFaceColor','y',...
            'MarkerSize',10); 
plot(c_star,Mc(c_star),'ko',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
legend('Using c','Using C*','quadrature points c','quadrature points C*')
hold off