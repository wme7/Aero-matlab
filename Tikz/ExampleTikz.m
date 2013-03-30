%% Example
% following ideas from http://almutei.wordpress.com
clear all; close all; clc;

% Do my Calculations
x=(0:0.25:7); sinx=sin(x); cosx=cos(x'); z = cosx*sinx;

% Plot my comparison of functions
figure(1)
hold on;
    h = plot(x,sinx,'--');
    k = plot(x,cosx,'-.');
hold off;
grid on; 
title('my comparison of sin(x) and cos(x)')
xlabel('x'); ylabel('y');
% convert to Tikz
matlab2tikz('myfigure.tex','height','4cm','width','3in')

% Plot surface
figure(2)
surf(z);
title('my surface z')
xlabel('x'); ylabel('y'); zlabel('z');
% convert to Tikz
matlab2tikz('myfigure2.tex','height','4cm','width','3in')