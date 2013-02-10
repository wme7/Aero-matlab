%% test ICs
clear all; clc; close all;

x = linspace(0,1,16);
uu = 1*(x<0.5)+0*(x>0.5);

nx = length(x);
dx = max(abs(x(1:end-1)-x(2:end)));
cfl = 0.6;
a = -1;
dt = dx*cfl/abs(a);
un = zeros(size(uu));

eps = 1E-25;
c1 = - 1/6;
c2 = 2/6;
c3 = 5/6;
c4 = - 7/6;
c5 = 11/6;

u = a*uu; % flux

%%
for i = 4:nx-3

B0p = 13/12*(u(i-2)-2*u(i-1)+u( i ))^2 + 1/4*(u(i-2)-4.*u(i-1)+3.*u(i))^2;
B1p = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-u(i+1))^2;
B2p = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(3.*u(i)-4.*u(i+1)+u(i+2))^2;

al0m = 1 / (10 * (eps + B0p))^2;
al1m = 6 / (10 * (eps + B1p))^2;
al2m = 3 / (10 * (eps + B2p))^2;

w0p = al0m / (al0m+al1m+al2m);
w1p = al1m / (al0m+al1m+al2m);
w2p = al2m / (al0m+al1m+al2m);

up(i) = w0p*(c1*u(i-2)+c3*u(i-1)+c2*u( i ))+...
        w1p*(c2*u(i-1)+c3*u( i )+c1*u(i+1))+...
        w2p*(c5*u( i )+c4*u(i+1)+c2*u(i+2));
end
           
h = up;

for j = 3:nx-4
    uun(j) = uu(j) - dt/dx*(h(j+1) - h(j));
end
