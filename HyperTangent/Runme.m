%% Playing with HyperTangent
% example on how to use the hyperbolic cut function in hcutfunc.m
clear all; close all; clc;

%% Define domain and ranges for the buffer region
x = 0:0.01:1; 
a = 0.45; b = 0.55;

% call hyperbolic cut function,
h = hcutfunc(x,a,b);

% Plot h in a figure
figure(1)
% plot shade
    x_shade = linspace(a,b,10);
    y1 = zeros(size(x_shade));
    y2 = 1.2*ones(size(x_shade));
    shade = shadedplot(x_shade,y1,y2,[0.7 0.7 1]);
    % [1 0.7 0.7]: red area
    % [0.7 0.7 1]: blue area
% plot h
hold on
    plot(x,h,'-.r');
    axis([0,1,0,1.2])
hold off
