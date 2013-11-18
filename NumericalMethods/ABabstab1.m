
% plot AB stability regions
theta = 0:0.001:2*pi;
z = exp(i*theta);

% AB2
% compute coefficients of integral of interpolant
A = [[ 1 -1];[1 0]];
cAB2 = [1 1/2]/A;

% compute margin of stability
nuAB2 = (z.^2-z)./(cAB2(2)*z+cAB2(1));

% AB3
% compute coefficients of integral of interpolant
A = [[1 -2 4];[1 -1 1];[1 0 0]];
cAB3 = [1 1/2 1/3]/A;

% compute margin of stability
nuAB3 = (z.^3-z.^2)./(cAB3(3)*z.^2 + cAB3(2)*z + cAB3(1));

% plot regions
plot(real(nuAB2),imag(nuAB2), 'r-');
hold on;
plot(real(nuAB3),imag(nuAB3), 'b-');

% plot axis
plot([0 0], [-10 10], 'k-');
plot([-10 10], [0 0], 'k-');
axis([-1.5 .5 -1 1]);
hold off;
legend('AB2', 'AB3')