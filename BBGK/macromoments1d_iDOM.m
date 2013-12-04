function [n, nux, E] = macromoments1d_iDOM(k,w,f,v,v_star)
% Interpolatory Discrete Ordinate Method for Computing integration moments

%% Inputs DATA
x = v; y = log(f);
[~,nx] = size(f);
[nc,~] = size(v_star);

% Uncomment to filter 'inf' values and identify usable data for interpolation
%y(isinf(y)) = 0; y_id = find(y);
%x = x(y_id); y = y(y_id); 
%[nv,nx] = size(y);

%% fitting: Least squares method
% Vandermonde matrix for a polynomial k=2 order,
V = [ones(size(x(:,1))) x(:,1) x(:,1).^2]; Vt = V';
a  = zeros(3,nx); for i = 1:nx; a(:,i) = ((Vt*V)^-1)*Vt*y(:,i); end;

%% Recover information from interpolation coefs,
%A2 = a(3); A1 = a(2); A0 = a(1);
sigma = sqrt(-1./(2*a(3,:)));
mu = a(2,:).*sigma.^2;
A = exp(a(1,:)+mu.^2./(2*sigma.^2));

%% Apply DOM
[sigma,mu,A] = apply_DOM(sigma,mu,A,nc);

%% Evaluating Interpolated funcion
% f_inter = A.*exp(-(x).^2);
f_inter = A.*exp(-v_star.^2);

%% Compute Integration Moments
% Following C.T. Hsu, et al.(2013)
c = mu + v_star.*sqrt((2*sigma.^2));
J = sqrt(2*sigma.^2); % Jacobian

% Using Quadrature rules to integrate for:
n   = k*sum(J.* w .* f_inter);    % Density
nux = k*sum(J.* c .* w .* f_inter);   % Density * velocity x
E   = k*sum(J.* 0.5.*( c.^2 ).* w .* f_inter);   % Total Energy Density
