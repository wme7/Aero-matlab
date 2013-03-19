clear all; close all; clc;
%Statistic
theta = 0; % MB / classical gas

% IC
rho = 1;
ux = 0.75;
p = 1;

% Thermodynamics tells us:
ne = p; % internal energy
%ne = 3/2*p; % internal energy

% Macroscopic props
%t = 4/3*ne./rho;
%z = rho./sqrt(pi.*t);
%p = 2/3*ne;
t = 4*ne./rho;
z = rho./sqrt(pi.*t);
p = ne;

% Gauss points and transformation
a = sqrt(t);
[c_star,w] = GaussHermite(3);
v = a*c_star + ux;  % transformation: v = a*(C*) + ux
w = w.*exp(c_star.^2);
J = a;

% Equilibrium
f = 1./(exp( (v-ux).^2 ./ t)./z + theta);

% Figure
plot(c_star,f,'-bo');

% Moments
rho_new = sum(w.*J.*f)
ux_new = sum(w.*J.*v.*f)./rho_new
ne_new = 1/2*sum(w.*J.*((v-ux).^2).*f)
p_new = ne_new
t_new = 4*ne./rho_new