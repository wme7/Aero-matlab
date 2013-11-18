
% plot real and imaginary parts of various finite difference derivative operators
% input: c = scalar multiplier 
% input: M = number of points on interval (excluding one end point)
% input: dx = uniform grid point separation          

function plotmultiplies(c, M, dx)

m=0:(M-1);

thetam = 2*pi*m/M;

deltaR = c*(exp(i*thetam)-1)/dx;
deltaL = c*(1-exp(-i*thetam))/dx;
deltaC = i*c*sin(thetam)/dx;
deltaC4 = i*c*(8*sin(thetam)-sin(2*thetam))/(6*dx);
deltaC6 = i*c*sin(thetam).*(1+(2/3)*sin(thetam/2).^2 +(8/15)*sin(thetam/2).^4)/dx;

% plot real part of eigenvalues
subplot(1,2,1); 
plot(thetam, real(deltaR), 'r-d');
hold on; 
plot(thetam, real(deltaL), 'b-*');
plot(thetam, real(deltaC), 'g-s');
plot(thetam, real(deltaC4), 'k-o');
hold off;
xlabel('\theta_m = 2\pim/M', 'FontSize', 18);
title('Real part of ODE multipliers', 'FontSize', 18);
legend('\delta_+', '\delta_-', '\delta_0', '\delta_0(1-(1/6)dx^2\delta^2)');

% plot imaginary part of eigenvalues 
subplot(1,2,2); 
plot(thetam, imag(deltaR), 'r-d');
hold on; 
plot(thetam, imag(deltaL), 'b-*');
plot(thetam, imag(deltaC), 'g-s');
plot(thetam, imag(deltaC4), 'k-o');
plot(thetam, imag(deltaC6), 'm-h');
plot([0 pi], [0 pi], 'k-', 'LineWidth', 2);      % first part of correct phase relationship
plot([pi 2*pi], [-pi 0], 'k-', 'LineWidth', 2); % second part

hold off;
xlabel('\theta_m = 2\pim/M', 'FontSize', 18);
title('Imaginary part of ODE multipliers', 'FontSize', 18);
legend('\delta_+', '\delta_-', '\delta_0', '\delta_0(1-(1/6)dx^2\delta^2)', '6th order central difference', 'Exact');




