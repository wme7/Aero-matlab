clear all; clc; close all;

% Parameteres,
u = 2.816; theta = 0; z = 0.515; t = 0.45;

% Initial Range for c,
c = -20:2:20;
c2 = -20:0.1:20; % for ploting 

% Generate sample data
f1 = 1./((1/z).*exp((c-u).^2/t)+theta);
f2 = 1./((1/z).*exp((c2-u).^2/t)+theta); % just for ploting

%% do x = c & y = log(f)
x = c;
y = log(f1);

% Filter inf values and identify usable data
y(isinf(y)) = 0; y_id = find(y);

% data for interpolation
x = x(y_id);
y = y(y_id);

%% Fitting
p = polyfit(x,y,2);
f3 = polyval(p,x);

% Data from coeficients
A2 = p(1); A1 = p(2); A0 = p(3);
sigma = sqrt(-1/(2*A2));
mu = A1*sigma^2;
A = exp(A0+mu^2/(2*sigma^2));

%% Compute GH information
[v_star,w] = GaussHermite(3); % for integrating range: -inf to inf
w = w.*exp(v_star.^2);

v = mu + v_star*sqrt((2*sigma^2));
%f4 = A*exp(-(v-mu).^2/(2*sigma^2));
f5 = A*exp(-v_star.^2);

rho = sum(sqrt(2*sigma^2)*w.*f5);

%% Visualization
subplot(1,2,1); grid on; hold on; 
                         plot(c2,f2,'--k','LineWidth',2); 
                         plot(v_star,f5,'og','LineWidth',2); 
                         hold off;
subplot(1,2,2); grid on; hold on; 
                         plot(x,y,'-*r','LineWidth',2); 
                         plot(x,f3,'-ob','LineWidth',2);
                         hold off;

% Using DOM for comparison:
g = @(x) 1./((1/z).*exp((x-u).^2/t)+theta);
tol = 1E-12;
rho_quad = quad(g,-20,20,tol);