%% Lagrange 
% Lagrage Interpolation Polinomial
clc;
clear;
close all;

%% Input Data
% Vector X
%x=[1 2 4 8];
%x=[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1];
x=[1978 1980 1982 1984 1986 1988 1990 1992];
%x=[1978 1980 1986 1988 1990 1992];
% Vector Y
%y=[1 3 7 11];
%y=[0.038 0.058 0.1 0.2 0.5 1 0.5 0.2 0.1 0.058 0.038];
y=[12 12.7 13 15.2 18.2 19.8 24.1 28.1];
%y=[12 12.7 18.2 19.8 24.1 28.1];
if length(x)==length(y)
    plot(x,y);
else
    fprintf('vectors must have the same length!')
end

%% Objetive value

%xi=7;
%xi=0.7;
%xi=1982;
%xi=1984;
xi=1994;

%% Lagrange Formulation
l=length(x);
lj=ones(1,l);
for j=1:1:l
    for n=1:1:l
        if j~=n
            lj(j)=lj(j)*(xi-x(n))/(x(j)-x(n));
        end
    end
end
lj;
yi=dot(lj,y)

%% Ploting Lagrange Result
% I wish to use 100 points for ploting the Lagrange polinomial.
%xpi=x(1):0.01:x(l);
xpi=x(1):0.01:xi;
lxpi=length(xpi);
for k=1:1:lxpi
    lj2=ones(1,l);
    for j=1:1:l
    for n=1:1:l
        if j~=n
            lj2(j)=lj2(j)*(xpi(k)-x(n))/(x(j)-x(n));
        end
    end
    end
    ypi(k)=dot(lj2,y);
end

%% Plot Comparison Graph
    plot(x,y,xpi,ypi);
