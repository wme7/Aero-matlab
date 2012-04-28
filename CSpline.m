%% Cubic Spline
% Cubic spline Interpolation 
clc;
clear;
close all;

%% Input Data
% Vector X
%x=[1 2 4 8];
%x=[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1];
%x=[1978 1980 1982 1984 1986 1988 1990 1992];
x=[1978 1980 1986 1988 1990 1992];
% Vector Y
%y=[1 3 7 11];
%y=[0.038 0.058 0.1 0.2 0.5 1 0.5 0.2 0.1 0.058 0.038];
%y=[12 12.7 13 15.2 18.2 19.8 24.1 28.1];
y=[12 12.7 18.2 19.8 24.1 28.1];
n=length(x);    %Number of elements to interpolate
if length(x)==length(y)
    plot(x,y);
else
    fprintf('vectors must have the same length!')
end

%% Objetive value

%xi=7;
%xi=0.7;
%xi=1982;
xi=1984;
%xi=1994;

%% Where is xi?
% find the number of the spicewise equation where yi would be computed
i=1;
no=0; %counter of equation
while i>=1
    if xi>x(n)
        no=n-1;
        i=0;
    elseif xi-x(i)>0
        no=no+1;
        i=i+1;
    else
        i=0;
    end
end

%% Forward Delta & Backward Delta
for i=1:1:n-1
    Deltaf(i)=x(i+1)-x(i);   %Forward Delta values
end
for i=2:1:n
    Deltab(i-1)=x(i)-x(i-1);  %Backward Delta values
end

%% Cubic Spline Formulation
% for A*g=f
% formulate vector f
for i=2:1:n-1
    f(i-1,1)=(y(i+1)-y(i))/Deltaf(i)-(y(i)-y(i-1))/Deltab(i);
end
% formulate Matrix A
A=zeros(n-2,n-2);   % matrix is always n-2, where n is the number of points to interpolate
A(1,1)=(Deltaf(1)+Deltab(1))/3; % for the first row elements
A(1,2)=Deltab(1)/6;
for i=2:1:n-3                    
    A(i,i-1)=Deltaf(i)/6;
    A(i,i)=(Deltaf(i)+Deltab(i))/3;
    A(i,i+1)=Deltab(i)/6;
end
A(n-2,n-3)=Deltaf(n-2)/6;       % for the last row elements
A(n-2,n-2)=(Deltaf(n-2)+Deltab(n-2))/3;
% formulate vector g
g=zeros(n,1);   %Column of unknown values (displacements)
go=inv(A)*f;
for i=1:1:n-2
    g(i+1)=go(i);
end

%% Boundary conditions
g(1)=0;         %free run-out condition 
g(n)=0;         %free run-out condition
%g(1)=g(2);      %Parabolic run-out condition
%g(n)=g(n-1);    %Parabolic run-out condition
%g(1)=g(n-1);    %Periodic run-out condition 
%g(n)=g(2);      %Periodic run-out condition

%% Interpolate yi
% we use the nearest equation gi(x) to compute yi:
i=no;
a=g(i)/6*((x(i+1)-xi)^3/Deltaf(i)-Deltaf(i)*(x(i+1)-xi));
b=g(i+1)/6*((xi-x(i))^3/Deltaf(i)-Deltaf(i)*(xi-x(i)));
c=y(i)/Deltaf(i)*(x(i+1)-xi)+y(i+1)/Deltaf(i)*(xi-x(i));
yi=a+b+c

xis=1978:0.1:1992;
yy=spline(x,y,xis);
h=plot(x,y,xis,yy);
