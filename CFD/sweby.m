%% Plot sweby diagrams

% number of points
N = 1000;
% compute flux limiting functions for a range of methods
theta = linspace(-1, 5, N);

upwind      = zeros(N,1);
laxwendroff = ones(N,1);
beamwarming = theta;
fromm       = 0.5*(1+theta);

minmo    = minmod(1, theta);
superbee = max(0, max(min(1,2*theta), min(2,theta)));
MC       = max(0, min(min( (1+theta)/2, 2), 2*theta));
vanleer  = (theta+abs(theta))./(1+abs(theta));

% plot TVD stability region 
% vertices of polygonal bounding region (truncated at x=10) 
figure(1)
x = [0 10 10 1];
y = [0   0   2 2];
fill(x,y,[.9 .9 .9]); % fill region with grey

hold on;
plot(theta, upwind,         'r-', 'LineWidth', 1.5)
plot(theta, laxwendroff,	'b-', 'LineWidth', 1.5)
plot(theta, beamwarming,    'g-', 'LineWidth', 1.5)
plot(theta, fromm,          'm', 'LineWidth', 1.5)
hold off;

axis equal
axis([-1 4 -.5 2.5])
legend('TVD region', 'Upwind', 'Lax-Wendroff', 'Beam-Warming', 'Fromm');
xlabel('\theta', 'FontSize', 18);
ylabel('\phi', 'FontSize', 18);

% plot Sweby modified TVD stability region 
% vertices of polygonal bounding region (truncated at x=10) 
figure(2)
x = [0 1 10 10 2  1 0.5 0];
y = [0 1  1  2 2  1   1 0];
fill(x,y,[.9 .9 .9]); % fill Sweby TVD region with grey

hold on;
plot(theta, minmo,      'r-', 'LineWidth', 1.5);
plot(theta, superbee,	'g-', 'LineWidth', 1.5);
plot(theta, MC,         'b-', 'LineWidth', 1.5);
plot(theta, vanleer,	'm-', 'LineWidth', 1.5);
hold off;

axis equal
axis([-1 4 -.5 2.5])
legend('Sweby TVD region', 'minmod', 'superbee', 'MC', 'van Leer');
xlabel('\theta', 'FontSize', 18);
ylabel('\phi', 'FontSize', 18);
