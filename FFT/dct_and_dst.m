% Using dcf and dst 
clear all; close all; clc;

f = @(x) x.^2/pi^2;

N = 16; j = 0:N;
x = (2*pi/N)*j;

y = f(x)';

yct = dct(y);
yst = dst(y);

figure(1)
hold on; grid on;
plot(abs(yst),'-or');
plot(abs(yct),'--sg');
legend('dst','dct')
hold off;