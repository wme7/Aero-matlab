clear all
% explicit euler
% m = 1; %adjust value of m
% y(1) = 1;%input your initial condition
% dt = 0.2; %adjust your step size
% T = 0:dt:15; %set up your time domain, here I have [0,5]
% for i = 2:length(T) %construct a 'for' loop
% y(i) = y(i-1) + m*dt*(-3*(i-2)/(i-1))*y(i-1)+m*dt*(2*(1+(i-2)^3))*exp(-(i-2));
% end
% plot(T,y)
% 
% Implicit euler
% n = 1; %adjust value of m
% w(1) = 1;%input your initial condition
% dtt = 0.2; %adjust your step size
% T1 = 0:dtt:15; %set up your time domain, here I have [0,5]
% for j = 2:length(T1) %construct a 'for' loop
%      w(j) = (w(j-1)+ dtt*(2*(j^3))*exp(-j+1))/(1+3*dtt*(j-1)/j) ;
% end
% plot(T1,w)
% 
% 
% Crank-Nicolson
m = 1; %adjust value of m
y(1) = 1;%input your initial condition
dt = 0.2; %adjust your step size
T = 0:dt:15; %set up your time domain, here I have [0,5]
for i = 2:length(T) %construct a 'for' loop
    y(i)=(dt*i^3*exp(-i+1)+dt*(i-1)^3*exp(-i+2)+...
           y(i-1)*(1-0.5*dt*(-3*i-6/(i-1))))/(1+0.5*dt*(3*i+3/i));
end
plot(T,y)





