% TEST
clear
%load('PDF'); load('PDF_eq'); load('v_DOM');
load('veq'); load('feq');

%% Inputs
y = log(f); x = v;
[nv,nx] = size(f);

%% fitting: Least squares method
% Vandermonde matrix for a polynomial order k=2 
V = [ones(size(x(:,1))) x(:,1) x(:,1).^2]; Vt = V';
a  = zeros(3,nx); for i = 1:nx; a(:,i) = ((Vt*V)^-1)*Vt*y(:,i); end;
b = [a(3,:);a(2,:);a(1,:)]; % order for polyval
% a = zeros(6,nx);
% for i = 1:nx; a(:,i) = polyfit(x(:,i),y(:,i),5); end;
% b = zeros(6,nx);
% for i = 1:nx; b(:,i) = polyfit(x(:,i),f(:,i),5); end;

%% Test Polynomial interpolation
y2 = zeros(nv,nx); for i=1:nx; y2(:,i) = polyval(b(:,i),x(:,i)); end;
%f_fit = zeros(nv,nx); for i=1:nx; f_fit(:,i) = polyval(b(:,i),x(:,i)); end;

%% Plot y1 and y2
subplot(2,2,1); surf(y);
subplot(2,2,2); surf(y2);
subplot(2,2,3); surf(exp(y));

%% Recover information from interpolation
A2 = a(3); A1 = a(2); A0 = a(1);
sigma = sqrt(-1./(2*a(3,:)));
mu = a(2,:).*sigma.^2;
A = exp(a(1,:)+ mu.^2./(2*sigma.^2));

% Apply DOM
[sigma,mu,A] = apply_DOM(sigma,mu,A,60);

%% Test proposal
f_next = A.*exp(-(x-mu).^2./(2*sigma.^2));
subplot(2,2,4); surf(f_next);
figure; surf(abs(f-f_next)) %compare output with original input
