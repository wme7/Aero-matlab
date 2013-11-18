% Adams-Bashforth

AB1 = [1];

AB2 = [3/2;-1/2];

AB3 = [23/12;-16/12;5/12];

AB4 = [55/24;-59/24;37/24;-9/24];


% AB Time stepping scheme
%
% y^{n+1} = y^{n} + mu*\sum_{m=1}^{m=r} ABr(m)*f^{n-m+1}
% 
% for a set of points on the unit circle we look for
% the required mu to keep the point fixed:
%
% y^n = z = e^{i*\theta}

theta = (0:0.01:2*pi)';
z = exp(sqrt(-1)*theta);

%
% AB1: mu = (z-1)/(1)
%

mu1 = (z-1);

%
% AB2: mu = 2*(z^2-z)/(3*z - 1)
%

mu2 = 2*(z.^2 - z)./(3*z-1);

%
% AB3: mu = 12*(z^3-z^2)/(23*z^2 - 16*z   + 5)
%

mu3 =  12*(z.^3-z.^2)./(23*z.^2 - 16*z   + 5);

%
% AB4: mu = 24*(z^4-z^3)/(55*z^3 - 59*z^2 + 37*z - 9)
%

mu4 = 24*(z.^4-z.^3)./(55*z.^3 - 59*z.^2 + 37*z - 9);

% plot stability limits

plot(real(mu1), imag(mu1), 'r-');
hold on; 
plot(real(mu2), imag(mu2), 'g-'); 
plot(real(mu3), imag(mu3), 'b-'); 
plot(real(mu4), imag(mu4), 'k-'); 
hold off;
legend('AB1', 'AB2', 'AB3', 'AB4');
xlabel('re(\mu)');
ylabel('im(\mu)');
axis equal;