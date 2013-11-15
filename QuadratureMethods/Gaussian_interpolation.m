%% Gaussian interpolation
clear all; %clc;

%initial parameters
z = 0.3;
theta = 1;

% Reference function
f = @(z,c,theta) 1./((1./z).*exp(c.^2) + theta);

% Compute reference data
x = -4:0.1:4;
fx = f(z,x,theta);
hold on
plot(x,fx,'--k');

% Testing values
xi = [-1.22474 0 1.22474]'; % GH abscissas
fi = f(z,xi,theta);
plot(xi,fi,'ob');

% fitting
ylog=log(fi);
xlog=xi;
p=polyfit(xlog,ylog,2);
A2=p(1); A1=p(2); A0=p(3);

sigma=sqrt(-1/(2*A2));
mu=A1*sigma^2;
A=exp(A0+mu^2/(2*sigma^2));
y=A*exp(-(x-mu).^2/(2*sigma^2));

plot(x,y,'or');
hold off

%% Compute derivate at xi points
dydxi = A*(mu-xi)/sigma^2.*exp(-(xi-mu).^2/(2*sigma^2));

%plot(x,dydxi(1)*x)
%plot(x,dydxi(2)*x)
%plot(x,dydxi(3)*x)