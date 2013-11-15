% Test 
clear all; clc; close all;

%Parameters
c = -7:7; z = 0.3; theta = 0; u = 1.2; t = 1.5;
g = @(c) 1./((1/z).*exp((c-u).^2/t) + theta); % Reference
f = g(c)';

% Reference integration
rho = quad(g,-7,7,1E-12);

nc = length(c);
% Legendre Vandermonde
[x,w,V] = GaussLegendre(nc);

% Transfrom to Legerdre Coefs,
ft = V\f;

f_sum = zeros(c,nc);
for j = x
    for i = 0:nc-1
        f_sum(j+8,i+1) = ft * legendreP(j,i);
    end
end
f_inter = sum(f_sum);

