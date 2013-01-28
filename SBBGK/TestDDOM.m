%% DDOM test
% This subroutine was abuilt to test a dynamic integration method

% IC's
theta = 0; z = 0.3; u = 0.0; t = 0.1;

% Create functions
M = @(c) 1./((1/z)*exp((c-u).^2/t) + theta); 
Mc = @(c_star) 1./((1/z)*exp(c_star.^2) + theta);

% Weight and C* values
nv = 3 %#ok<NOPTS>
[c_star,wc] = GaussHermite(nv);
[v,wv] = GaussHermite(nv);

% Transformation rule of c and C*
c = sqrt(t) * c_star + u; 

% Evaluate
feq  = M(v);
feqc = Mc(c_star);
if feq == feqc; fprintf('feq & feqc agree \n\n'); else fprintf('feq & feqc disgree \n\n'); end;

%% 1D  CASE
J = sqrt(t);

% GH quadrature function weight:
wv = wv.*exp(v.^2); 
wc = wc.*exp(c_star.^2);
fprintf('\n');

% Apply GH quadrature rule:
tic; n1 = sum(0.5*(v-u).^2'.*wv'.*feq'); toc; 
tic; n2 = sum(0.5*(c-u).^2'.*wc'.*feqc'.*J); toc; 
fprintf('n1: %1.12f \n\n',n1);
fprintf('n2: %1.12f \n\n',n2); 
if n1 == n2; fprintf('n1 & n2 agree \n\n'); else fprintf('n1 & n2 disgree \n\n'); end;

% reference solution
tol = 1E-12;
tic; n_ref = quad(M,-7,7,tol); toc; 
fprintf('n_ref: %1.12f \n\n',n_ref);

% plot
hold on
x = -4:0.1:4; grid on;
plot(x,M(x)); plot(x,Mc(x));
hold off