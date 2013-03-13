%% Dynamic Discrete Ordinate Method (DDOM)
clear all; close all; clc;

% Example function of the Equilibrium Distribution Function
origin = 0;
switch origin
    case 0 % centered
        z = 0.30; disp = 0.0; % normalization constant or fugacity
        theta = -1; F0_BE = @(c) 1./((1./z).*exp(c.^2) + theta);
        theta =  0; F0_MB = @(c) 1./((1./z).*exp(c.^2) + theta);
        theta =  1; F0_FD = @(c) 1./((1./z).*exp(c.^2) + theta);
    case 1 % displaced
        z = 0.10; disp = 2.0; t = 3.2;
        theta = -1; F0_BE = @(c) 1./((1./z).*exp((c-disp).^2) + theta);
        theta =  0; F0_MB = @(c) 1./((1./z).*exp((c-disp).^2) + theta);
        theta =  1; F0_FD = @(c) 1./((1./z).*exp((c-disp).^2) + theta);
    otherwise
        error('invalid case')
end

% Discretize velocity domain
a = -7; b = 7; c = a:0.1:b;

% Compute velocity distributions for every statistic
f_BE = F0_BE(c);
f_MB = F0_MB(c);
f_FD = F0_FD(c);

% Integrate for density values using Matlab's ode45
tol = 10E-10; % Tolerance
tic; rho_BE = quad(F0_BE,a,b,tol); toc; fprintf('rho_BE: %1.12f\n\n',rho_BE) 
tic; rho_MB = quad(F0_MB,a,b,tol); toc; fprintf('rho_MB: %1.12f\n\n',rho_MB)
tic; rho_FD = quad(F0_FD,a,b,tol); toc; fprintf('rho_FD: %1.12f\n\n',rho_FD)

% Gauss Hermite Weights and Abscissas
nv = 7; [v,w] = GaussHermite(nv); w = w.*exp(v.^2);

% Integrate for density values using Gauss Hermite Quadrature
tic; rho2_BE = sum(F0_BE(v).*w); toc; fprintf('rho_BE: %1.12f\n\n',rho2_BE)
tic; rho2_MB = sum(F0_MB(v).*w); toc; fprintf('rho_MB: %1.12f\n\n',rho2_MB)
tic; rho2_FD = sum(F0_FD(v).*w); toc; fprintf('rho_FD: %1.12f\n\n',rho2_FD)

% Error
err_BE = abs(rho_BE-rho2_BE); fprintf('GH quad in BE error = %e\n\n',err_BE);
err_MB = abs(rho_MB-rho2_MB); fprintf('GH quad in MB error = %e\n\n',err_MB);
err_FD = abs(rho_FD-rho2_FD); fprintf('GH quad in FD error = %e\n\n',err_FD);

% Plot shapes of three statistics
range = [a,b,0,max(f_BE)]; 
subplot(2,2,1); plot(c,f_BE); title('BE'); 
    axis(range); xlabel('nv'); ylabel('f'); grid on;
subplot(2,2,2); plot(c,f_MB); title('MB'); 
    axis(range); xlabel('nv'); ylabel('f'); grid on;
subplot(2,2,3); plot(c,f_FD); title('FD'); 
    axis(range); xlabel('nv'); ylabel('f'); grid on;
subplot(2,2,4),axis('off'), title('GH Quadrature Error','fontsize',12)
    text(0.0,0.95,'Inputs for f^{eq}_{quantum}','fontsize',12)
    text(0.1,0.80,['z = ',num2str(z)],'fontsize',12)
    text(0.1,0.65,['a = ',num2str(disp)],'fontsize',12)
    text(0.0,0.50,'Output: Numerical Integration error','fontsize',12)
    text(0.1,0.35,['Error in BE = ',num2str(err_BE)],'fontsize',12)
    text(0.1,0.20,['Error in MB = ',num2str(err_MB)],'fontsize',12)
    text(0.1,0.05,['Error in FD = ',num2str(err_FD)],'fontsize',12)