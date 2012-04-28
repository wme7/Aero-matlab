clear all
close all
clc

x=0:0.01:pi;
w=0:0.02:2*pi;
y=w.*sin(x);

hold on
plot(x,y','red');
plot(x,w,'yellow');
title('example plot')
grid on
xlabel('x')
ylabel('sin(x)')
legend('my function','w')
hold off