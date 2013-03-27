%% Laboratory No.3, Part 2: Solving for the Tappered beam.
% Coded by Manuel Diaz, 2013.03.20
clear all; close all; clc;

% Load cases per number of elements
elements = [1 2 4 8 16];

% Solve and plot Displacement
figure(1)
hold on
for j = 1:5
   [D,R,S,x] = Driver1d(elements(j));
   markers = {':+r','--og','-.*b',':xm','--sk','-.dw'}; 
   plot(x,D,markers{j})
end
title('Displacement vs x'); ylabel('displacment (m)'); xlabel('x (m)');

% Compute exact solution
exact_D= 1/35000*log(x+0.5);
plot(x,exact_D,'-r')
legend('1 element','2 element','4 element','8 element','16 element','exact',2)
hold off

% Compute and plot stress
figure(2)
hold on
for j = 1:5
   [D,R,S,x] = Driver1d(elements(j));
   markers = {':+r','--og','-.*b',':xm','--sk','-.dw'}; 
   S = [S;S(end)]; % repeat final value to match number of nodes
   stairs(x,S,markers{j})
   clear S
end
title('Stress vs x'); ylabel('displacment (m)'); xlabel('x (m)');

% Compute exact solution
exact_S = 1/35000*70e9./(x+0.5);
plot(x,exact_S,'-r')
legend('1 element','2 element','4 element','8 element','16 element','exact',1)
hold off
