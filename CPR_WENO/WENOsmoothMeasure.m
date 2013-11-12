%% Just for testing
tic
clear all; clc;

K = 5;

xgrid = mesh1d([0,1],10,'LGL',K);

w = xgrid.weights';
x = xgrid.nodeCoordinates;
dx = xgrid.elementSize;
xi = xgrid.solutionPoints;

u = IC(x,3);

u_bar = w*u/2;

L = LagrangePolynomial(xi);

%% Trouble cells
% assume we now the troubled Elements
troubleE = 5:7;
toc % 0.06 sec -> I can accept this.

%% 
disp('Compute Beta coefs')
tic
Bcoef = zeros(K,K+1);
for s = 1:K 
    syms x;
    dpdx  = L.dnlagrangePolynomial(s);
    Bcoef(s,:) = double(int((dx)^(2*s-1)*(dpdx).^2,x,-1,1));
end

% Beta factors for every element
B = sum(Bcoef*u);
toc % 2 sec! -> must be prebuilt

%%
tic
gamma = [1e-6,0.999998,1e-6];
epsilon = 1e-6;
w_tilde = gamma./(epsilon+B(troubleE)).^2;
w0 = w_tilde(1)/sum(w_tilde);
w1 = w_tilde(2)/sum(w_tilde);
w2 = w_tilde(3)/sum(w_tilde);
toc % 0.001 sec

%%
disp('Exact WENO limiter implementation')
tic
u_new = zeros(size(u));
for j = troubleE
    P0 = L.lagrangePolynomial*u(:,j-1);
    P1 = L.lagrangePolynomial*u(:,j);
    P2 = L.lagrangePolynomial*u(:,j+1);
    P0_bar = u_bar(j-1);
    P1_bar = u_bar(j);
    P2_bar = u_bar(j+1);
    P0_tilde = P0 - P0_bar + P1_bar;
    P2_tilde = P2 - P2_bar + P1_bar;
    u_new(:,j) = subs(w0*P0_tilde,xi+2) + subs(w1*P1,xi) + subs(w2*P2_tilde,xi-2);
end
toc
% 3 sec! -> no f***ing way!!!

%% build polynomial function 
% for the a Lagrange Polynomial in the standard element
disp('Approximate WENO limiter implementation')
tic
switch K
    case 3
        phi0 = LGL_K3(xi+2);
        phi1 = LGL_K3(xi);
        phi2 = LGL_K3(xi-2);
    case 4
        phi0 = LGL_K4(xi+2);
        phi1 = LGL_K4(xi);
        phi2 = LGL_K4(xi-2);
    case 5
        phi0 = LGL_K5(xi+2);
        phi1 = LGL_K5(xi);
        phi2 = LGL_K5(xi-2);
end
toc

tic
u_new2 = zeros(size(u));
for j = troubleE
    P0 = phi0*u(:,j-1);
    P1 = phi1*u(:,j);
    P2 = phi2*u(:,j+1);
    P0_bar = u_bar(j-1);
    P1_bar = u_bar(j);
    P2_bar = u_bar(j+1);
    P0_tilde = P0 - P0_bar + P1_bar;
    P2_tilde = P2 - P2_bar + P1_bar;
    u_new2(:,j) = w0*P0_tilde + w1*P1 + w2*P2_tilde;
end
toc
%fair result

%%
plotrange = [0,1,0.9,2.1]; 
subplot(2,2,1)
plot(xgrid.nodeCoordinates,u); axis(plotrange); grid on;
title('Initial Condition'); xlabel('x'); ylabel('u(x,t)');
subplot(2,2,2)
stairs(xgrid.elementCenter,u_bar); axis(plotrange); grid on;
title('Cell Averages in E_j'); xlabel('x'); ylabel('u(x,t)');
subplot(2,2,3)
plot(xgrid.nodeCoordinates,u_new); axis(plotrange); grid on;
title('Exact Limiter Implementation'); xlabel('x'); ylabel('u(x,t)');
subplot(2,2,4)
plot(xgrid.nodeCoordinates,u_new2); axis(plotrange); grid on;
title('Approximate Limiter Implementation'); xlabel('x'); ylabel('u(x,t)');
