function [n, j_x, j_y, E] = macromoments2d(k,w1,w2,f,v1,v2)
%% Compute Macroscopic Moments
% Precompute this product just to speed up things
w12 = w1 .* w2;

% Using Quadrature rules to integrate for:
n(:,:)   = k*k*sum(sum(w12 .* f));    % Number density
j_x(:,:) = k*k*sum(sum(v1 .* w12 .* f));   % Macrospic moment in x
j_y(:,:) = k*k*sum(sum(v2 .* w12 .* f));   % Macrospic moment in y
E(:,:)   = k*k*sum(sum(1/2*( v1.^2 + v2.^2 ).* w12 .* f)); % Energy Density