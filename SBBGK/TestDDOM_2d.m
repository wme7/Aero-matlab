%% DDOM test
% This subroutine was built to test a dynamic integration method in 2d.
clear all; clc; close all;

% IC's
%theta = 0; z = 0.3000; ux = 1.00; uy = 1.00; t = 0.1; % random
%theta = 0; z = 0.0394; ux = 0.00; uy = 0.00; t = 3.2; % Sod's right
theta = 0; z = 0.2821; ux = 0.75; uy = 0.75; t = 4.0; % Sod's left

% Create functions
M = @(cx,cy) 1./((1/z)*(exp( ((cx-ux).^2 + (cy-uy).^2 ) ./t) + theta));
Mc = @(cx_star,cy_star) 1./((1/z)*exp(cx_star.^2 + cy_star.^2) + theta);

% Weight and C* values
nv = input('Number of quadrature desided quadrature points: ');
%nv = 3; 
fprintf('quadrature points used: %1.0f x %1.0f \n\n',nv,nv);

[c_star,wc] = GaussHermite(nv);
wc = wc.*exp(c_star.^2);
[cx_star,cy_star] = meshgrid(c_star,c_star');
[wcx,wcy] = meshgrid(wc,wc');
    
[v,wv] = GaussHermite(nv);
wv = wv.*exp(v.^2);
[vx,vy] = meshgrid(v,v');
[wvx,wvy] = meshgrid(wv,wv');
    
%% 2D  CASE
% Transformation rule of c and C*
cx = sqrt(t) * cx_star + ux; 
cy = sqrt(t) * cy_star + uy; 

% Jacobian
J = sqrt(t); % J = a

% Evaluate
feq  = M(vx,vy) %#ok<NOPTS>
feqc = Mc(cx_star,cy_star) %#ok<NOPTS>
fprintf('\n');

%% Apply GH quadrature rule:
fprintf('***************************************************************\n')
fprintf('int(f_eq) \n\n')
tic; n1 = sum(sum(wvx.*wvy.*feq)); toc; 
tic; n2 = sum(sum(wcx.*wcy.*feqc.*J.*J)); toc; 
fprintf('\n int(f_eq): %1.12f \n',n1);
fprintf(' int(f_eqc*J): %1.12f \n\n',n2); 

fprintf('***************************************************************\n')
fprintf('int(cx*f_eq) \n\n')
tic; n3 = sum(sum(vx.*wvx.*wvy.*feq)); toc; 
tic; n4 = sum(sum(cx.*wcx.*wcy.*feqc.*J.*J)); toc; 
fprintf('\n int(cx*f_eq): %1.12f \n',n3);
fprintf(' int(cx*f_eqc*J): %1.12f \n\n',n4);

fprintf('***************************************************************\n')
fprintf('int(cy*f_eq) \n\n')
tic; n5 = sum(sum(vy.*wvx.*wvy.*feq)); toc; 
tic; n6 = sum(sum(cy.*wcx.*wcy.*feqc.*J.*J)); toc; 
fprintf('\n int(cy*f_eq): %1.12f \n',n5);
fprintf(' int(cy*f_eqc*J): %1.12f \n\n',n6);

fprintf('***************************************************************\n')
fprintf('int(1/2*((cx-ux)^2+(cy-uy)^2)*f_eq) \n\n')
tic; n7 = sum(sum(0.5*((vx-ux).^2+(vy-uy).^2).*wvx.*wvy.*feq)); toc; 
tic; n8 = sum(sum(0.5*((cx-ux).^2+(cy-uy).^2).*wcx.*wcy.*feqc.*J.*J)); toc; 
fprintf('\n int(1/2*((cx-ux)^2+(cy-uy)^2)*f_eq): %1.12f \n',n7);
fprintf(' int(1/2*((cx-ux)^2+(cy-uy)^2)*f_eqc*J): %1.12f \n\n',n8);

%% Reference solution
fprintf('***************************************************************\n')
fprintf('Using "quad2d" function \n\n')
%tol = 1E-12;
M1 = @(x,y) 1./((1./z).*(exp( ((x-ux).^2 + (y-uy).^2 ) ./t) + theta));
tic; n_ref = quad2d(M1,-7,7,-7,7); toc; 
fprintf('int(f_eq): %1.12f \n\n',n_ref);
M2 = @(x,y) x./((1./z).*(exp( ((x-ux).^2 + (y-uy).^2 ) ./t) + theta));
tic; n_ref = quad2d(M2,-7,7,-7,7); toc; 
fprintf('int(cx*f_eq): %1.12f \n\n',n_ref);
M3 = @(x,y) y./((1./z).*(exp( ((x-ux).^2 + (y-uy).^2 ) ./t) + theta));
tic; n_ref = quad2d(M3,-7,7,-7,7); toc; 
fprintf('int(cy*f_eq): %1.12f \n\n',n_ref);
M4 = @(x,y) (0.5*((x-ux).^2+(y-uy).^2))./((1./z).*(exp( ((x-ux).^2 + (y-uy).^2 ) ./t) + theta));
tic; n_ref = quad2d(M4,-7,7,-7,7); toc; 
fprintf('int(1/2*((cx-ux)^2+(cy-uy)^2)*f_eq): %1.12f \n\n',n_ref);
fprintf('***************************************************************\n')

%% Plot comparison solution
xx = linspace(-7,7,70);
yy = xx';
[y,x] = meshgrid(yy,xx);
subplot(1,2,1)
surf(M(x,y)); title('using discrete velocities cx and cy');
set(gca,'xDir','reverse');
   xlabel('vx - Velocity Space'); 
   ylabel('vy - Velocity Space');
   zlabel('f - Probability');
subplot(1,2,2)
surf(Mc(x,y)); title('using discrete velocities Cx* and Cy*');
set(gca,'xDir','reverse');
   xlabel('cx - Velocity Space'); 
   ylabel('cy - Velocity Space');
   zlabel('f - Probability');