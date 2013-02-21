% Using Chebyshev Polynomials
clear all; close all; clc;

% Define:
f = @(x) x.^4;
%f = @(x) 4*(x.^2-x.^4).*exp(-x./2);

x = (-1:0.1:1)'; % for -1 <= x <= 1

