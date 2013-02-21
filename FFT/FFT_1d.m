% Using FFT 
clear all; close all; clc;

f = @(x) cos(3*x);
%f = @(x) 1*(x>=0 & x<pi) + (-1)*(x>=pi & x<2*pi);

N = 16; j = 0:N-1;

x = (2*pi/N)*j;

y = f(x)';

yt = fft(y);

% Plot coeficients
figure(1)
grid on; plot(abs(yt),'or'); legend('fft')
