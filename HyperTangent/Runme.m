%% Playing with HyperTangent
% example on how to use the hyperbolic cut function in hcutfunc.m
clear all; close all; clc;

%% Define domain and ranges for the buffer region
x = 0:0.01:1; 
a = 0.4; b = 0.6;

%% Choose cut function model
cutmodel = 1; % {1} scaled tanh(x), {2} scaled cos(x), {3} linear model

%% Call cut function,
h = hcutfunc(x,a,b,cutmodel);

%% Plot h(x) in a figure
figure(1)
% plot shade
    x_shade = linspace(a,b,10);
    y1 = -0.2*ones(size(x_shade));
    y2 =  1.2*ones(size(x_shade));
    [shade,l1,l2] = shadedplot(x_shade,y1,y2,[0.7 0.7 1]);
    % [1 0.7 0.7]: red area
    % [0.7 0.7 1]: blue area
% plot h
hold on
    func = plot(x,h,'-.k','LineWidth',2);
    legend([shade(2),func],{'Buffer region','Cut function h(x)'},2);
    xlabel('x'); ylabel('h(x)');
    axis([0,1,-0.2,1.2]);
hold off
