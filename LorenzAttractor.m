%% Lorenz Attractor
% Solving a Non-Linear System of Differential Equations with several
% degrees of freedom.
%
% $\frac{\partial x}{\partial t} = \sigma (y-x)$ , 
% $\frac{\partial y}{\partial t} = rx-y-xz$ and  
% $\frac{\partial z}{\partial t} = xy-bz$
%
% where the values of $\sigma$ and $b$ are fixed leaving $r$ as the control
% parameter.
%
% For low values of $r$, the stable solution are stationary. When r exceeds
% 24.74, the trayectories in $x y z$ space become irregular orbits about two
% particular points.

%% Fixed Values
sigma = 10; b =8/3;

%% [Case 1] 
% Solve these equations for $r = 20$, starting from point P(1,1,1)

% Controling parameter
r = 20;
% Initial Conditions
x0 = 1; y0 = 1; z0 = 1;

%% Solve ODEs for [Case 1]

ivp = [x0,y0,z0,sigma,b,r];
%options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,u] = ode45(@LorenzAttractorfunc,[0 25], ivp);
x = u(:,1); y = u(:,2); z = u(:,3);

%% Make preatty figures for [Case 1]
% 1. plot the trajectory for 0<t<25 in the xy, yz and xz planes
figure
subplot(2,2,1)
plot(x,y);
title('xy trajectory');
xlabel('-x-');ylabel('-y-');

subplot(2,2,2)
plot(x,z);
title('xz trajectory');
xlabel('-x-');ylabel('-z-');

subplot(2,2,3)
plot(y,z);
title('yz trajectory');
xlabel('-y-');ylabel('-z-');

% 2. Plot x, y and z vs t and comments your plots
figure

subplot(3,1,1)
plot(t,x,'green'); 
title('x with r = 20'); 
ylabel('-x-'); xlabel('-t-');  
grid on

subplot(3,1,2)
plot(t,y,'blue'); 
title('y with r = 20'); 
ylabel('-y-'); xlabel('-t-');  
grid on

subplot(3,1,3)
plot(t,z,'red'); 
title('z with r = 20'); 
ylabel('-z-'); xlabel('-t-');  
grid on

% 3. for this case, plot the solution in xyz space, letting $z$ be your
% horizontal plane
figure; plot3(z,x,y);
xlabel('-z-'); ylabel('-x-'); zlabel('-y-');
title('x,y,z trajectory at r=20');
grid on

%% [Case 2] 
% Observe the change in the solution by repeating Case 1 for $r = 28$

%Controling parameter
r = 28;
% Initial Conditions
x0 = 1; y0 = 1; z0 = 1;

%% Solve ODEs for [Case 2]

ivp2 = [x0,y0,z0,sigma,b,r];
%options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t2,w] = ode45(@LorenzAttractorfunc,[0 25], ivp2);
x = w(:,1); y = w(:,2); z = w(:,3);

%% Make preatty figures for [Case 2]
% 1. plot the trajectory for 0<t<25 in the xy, yz and xz planes
figure
subplot(2,2,1)
plot(x,y);
title('xy trajectory');
xlabel('-x-');ylabel('-y-');

subplot(2,2,2)
plot(x,z);
title('xz trajectory');
xlabel('-x-');ylabel('-z-');

subplot(2,2,3)
plot(y,z);
title('yz trajectory');
xlabel('-y-');ylabel('-z-');

% 2. Plot x, y and z vs t and comments your plots
figure

subplot(3,1,1)
plot(t2,x,'green'); 
title('x with r = 28'); 
ylabel('-x-'); xlabel('-t-');  
grid on

subplot(3,1,2)
plot(t2,y,'blue'); 
title('y with r = 28'); 
ylabel('-y-'); xlabel('-t-');  
grid on

subplot(3,1,3)
plot(t2,z,'red'); 
title('z with r = 28'); 
ylabel('-z-'); xlabel('-t-');  
grid on

% 3. for this case, plot the solution in xyz space, letting $z$ be your
% horizontal plane
figure; plot3(z,x,y);
xlabel('-z-'); ylabel('-x-'); zlabel('-y-');
title('x,y,z trajectory at r=28');
grid on

%% [Case 3] 
% Observe the unpredictability at $r=28$ by overploting two solutions
% versus time starting from two initial nearby points:
% point 1: (6,6,6)
% point 2: (6,6.01,6)

%Controling parameter
r = 28;
% Initial Conditions
ivp3  = [6,6,6,sigma,b,r];
ivp4  = [6,6.01,6,sigma,b,r];

%% Solve ODEs for [Case 3]

options = odeset('RelTol',1e-14,'AbsTol',1e-14);
[time1,v] = ode45(@LorenzAttractorfunc,[0 25], ivp3, options);
x1 = v(:,1); y1 = v(:,2); z1 = v(:,3);

options = odeset('RelTol',1e-14,'AbsTol',1e-14);
[time2,w] = ode45(@LorenzAttractorfunc,[0 25], ivp4, options);
x2 = w(:,1); y2 = w(:,2); z2 = w(:,3);

%% Make preatty figures for [Case 3]

% 1. Plot x, y and z vs t and comments your plots
figure

subplot(3,1,1)
hold on
plot(time1,x1,'green'); plot(time2,x2,'blue');
title('x with r = 28'); 
legend('from P(6,6,6)','from P(6,6.01,6)');
xlabel('-t-'); 
grid on;
hold off

subplot(3,1,2)
hold on
plot(time1,y1,'green'); plot(time2,y2,'blue');
title('y with r = 28'); 
legend('from P(6,6,6)','from P(6,6.01,6)');
xlabel('-t-'); 
grid on;
hold off

subplot(3,1,3)
hold on
plot(time1,z1,'green'); plot(time2,z2,'blue');
title('z with r = 28'); 
legend('from P(6,6,6)','from P(6,6.01,6)');
xlabel('-t-'); 
grid on;
hold off

% 2. Compare both trajectories in xyz space:
figure
plot3(z1,x1,y1,'green',z2,x2,y2,'blue');
legend('from P(6,6,6)','from P(6,6.01,6)');
xlabel('-z-'); ylabel('-x-'); zlabel('-y-');
title('z,x,y trajectory at r=28');
grid on