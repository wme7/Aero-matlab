%% Example
% following ideas from http://almutei.wordpress.com
clear all; close all; clc;

% Do my Calculations
x=(0:0.2:2*pi); sinx=sin(x); cosx=cos(x'); z = cosx*sinx;

% Plot my comparison of functions
figure(1)
hold on;
    h = plot(x,sinx,'--r');
    k = plot(x,cosx,'-.g');
hold off;
grid on; % note that this comand is ignored by tikz
title('my comparison of sin(x) and cos(x)')
box on; xlabel('x'); ylabel('y'); 
% Convert to Tikz
matlab2tikz('myfigure.tex','height','4cm','width','3in')

% Plot surface
figure(2)
mesh(z);
title('my surface z')
xlabel('x'); ylabel('y'); zlabel('z');
% Convert to Tikz
matlab2tikz('myfigure2.tex','height','4cm','width','3in')

% Plot Contour
figure(3)
contour(z);
title('my surface z')
xlabel('x'); ylabel('y'); zlabel('z');
% Convert to Tikz
matlab2tikz('myfigure3.tex','height','4cm','width','3in')