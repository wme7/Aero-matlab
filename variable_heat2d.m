%% Heat Equation with variable coeficients
% Solution of the heat equation using point gauss method
%
% $$\frac{\partial}{\partial x}\left (k(x,y)\frac{\partial{T}}{\partial{x}} \right )+\frac{\partial}{\partial x}\left (k(x,y)\frac{\partial{T}}{\partial{y}} \right)=0$$
%
% where:
%
% $$k(x,y)=1+exp\left[\left( \frac{x-0.25-0.5y}{0.25} \right)^2 \right]$$
clear
clc
close all

%% Discretization of space
% Length and heigth of our 2D domain:
L = 1; % length in x direction
H = 1; % height in y direction

% our target is to evaluate the step h size, thus
h = 0.05; % we asume and equally spaced XY grid
hx = h;
hy = h; 

% number of points
n = round(L/hx);
m = round(H/hy);

% initalize Matriz T and k
T=zeros(n+1,m+1);
k=zeros(n+1,m+1);

% initalize Matriz k(x,y)
for i=1:n+1
    for j=1:m+1
        k(i,j)=1+exp(-((i/(n+1)-0.25-0.5*j/(m+1))/0.25)^2);
    end
end

%% Boundary Conditions
% Bottom (last row of the matrix)

T(1,:)=0;

% Top (first row of the matrix)
%T(1,:)=1;
for i=0:n
    % using function of x to prescribe the values.
    T(m+1,i+1)=1-cos(pi/(n)*i); 
end

% Left and Right
%Both conditions would be updated as the computation 

%% Method of Solution
%
T1=T;
for s=1:1000
    for i=2:n
        for j=2:m
            A(i,j)=k(i+1,j)+4*k(i,j)-k(i-1,j);
            B(i,j)=k(i,j+1)+4*k(i,j)-k(i,j-1);
            C(i,j)=k(i-1,j)+4*k(i,j)-k(i+1,j);
            D(i,j)=k(i,j-1)+4*k(i,j)-k(i,j+1);
            T1(i,j)=1/(4*k(i,j))*(0.25*A(i,j)*T1(i+1,j)+0.25*B(i,j)*T1(i,j+1)+0.25*C(i,j)*T1(i-1,j)+0.25*D(i,j)*T1(i,j-1));
            T1(2:m,1)=T1(2:m,2);
            T1(2:m,n+1)=T1(2:m,n);
        end
    end
end

%% Plot

%surf(T1)
contourf(T1)
colormap hot
colorbar('location','southoutside')

%% Method of Solution 2
%
T2=T;
for s=1:1000
    for i=2:n
        for j=2:m
            A(i,j)=k(i,j)+k(i+1,j);
            B(i,j)=k(i,j)+k(i,j+1);
            C(i,j)=k(i,j)+k(i-1,j);
            D(i,j)=k(i,j)+k(i,j-1);
            T2(i,j)=1/(k(i+1,j)+k(i-1,j)+4*k(i,j)+k(i,j+1)+k(i,j-1))*(A(i,j)*T2(i+1,j)+B(i,j)*T2(i,j+1)+C(i,j)*T2(i-1,j)+D(i,j)*T2(i,j-1));
            T2(2:m,1)=T2(2:m,2);
            T2(2:m,n+1)=T2(2:m,n);
        end
    end
end

%% Plot

%surf(T2)
contourf(T2)
colormap hot
colorbar('location','southoutside')
