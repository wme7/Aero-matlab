
function plotexp(mu)

time = 0:.001:10;

f = exp(mu*time);

plot(time, real(f), 'r-');
hold on;
plot(time, imag(f), 'b-');
plot(time, abs(f), 'k-');

nre = sprintf('re(exp(\mu t))');
nim = sprintf('im(exp(\mu t))');
nab = sprintf('|exp(\mu t)|');
legend('re(exp(\mut))', 'im(exp(\mut))', '|exp(\mut)|');

title(sprintf('mu = %g + i%g', real(mu), imag(mu)), 'FontSize', 18)
xlabel('t', 'FontSize', 18);
hold off;
rmax = max([real(f),abs(f),imag(f)]);
rmin = min([real(f),abs(f),imag(f)]);

rrange = rmax-rmin;
rmin = rmin-.1*rrange;
rmax = rmax+.1*rrange;

axis([-.1 10.1 rmin rmax])

