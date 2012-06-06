%% Homework 12 (b)
% Heat equation with variable coeficients using Gauss-Seidel Method:
%
% $\frac{\partial}{\partial x}\left (k(x,y)\frac{\partial{T}}{\partial{x}} \right )+\frac{\partial}{\partial x}\left (k(x,y)\frac{\partial{T}}{\partial{y}} \right)=0$
%
% where:
%
% $k(x,y)=1+exp\left[-\left( \frac{x-0.25-0.5y}{0.25} \right)^2 \right]$

%close all
clear all

%% Grid
h = 0.02;
lx = [0,1]; dx = h; x = lx(1):dx:lx(2); n = length(x);
ly = [0,1]; dy = h; y = ly(1):dy:ly(2); m = length(y);
T = zeros(m,n); T_next = zeros(m,n);
 
% Computing k(x,y) on the grid
k = zeros(m,n);
for j=1:m
    for i=1:n
        %k(j,i)=1+exp(-((x(i)-0.25-0.5*y(j))/0.25)^2);
        k(j,i)=1;
    end
end

%% Boundary Conditions
T(1,:) = 0; % Bottom (last row of the matrix)
T(m,:) = 1-cos(pi*x); % Top (first row of the matrix)
% Left and Right BCs will be updated on every iteration 

%% Method of Solution
%
A = zeros(m,n); B = zeros(m,n); C = zeros(m,n); D = zeros(m,n);
r = 0.00001;
iter = 0;
r_iter = 1;
while r_iter >= r
%for s=1:3000
    % we know the boundary values from the begining
    T_next(1,:) = 0;            % matrix bottom row, Dirichlet BC
    T_next(m,:) = 1-cos(pi*x);  % matrix top row, Dirichlet BC
    for i=2:m-1
        for j=2:n-1
            A(i,j)=k(i,j)+k(i+1,j);
            B(i,j)=k(i,j)+k(i,j+1);
            C(i,j)=k(i,j)+k(i-1,j);
            D(i,j)=k(i,j)+k(i,j-1);
            if j == 2 % Matrix left Column, Neumann BC
                T_next(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+3*k(i,j)+k(i,j+1))*(...
                    A(i,j)*T(i+1,j)+...
                    B(i,j)*T(i,j+1)+...
                    C(i,j)*T_next(i-1,j));
            elseif j == n-1 % Matrix right Column, Neumann BC
                T_next(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+3*k(i,j)+k(i,j-1))*(...
                    A(i,j)*T(i+1,j)+...
                    C(i,j)*T_next(i-1,j)+...
                    D(i,j)*T_next(i,j-1));
            else
                T_next(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+4*k(i,j)+k(i,j+1)+k(i,j-1))*(...
                    A(i,j)*T(i+1,j)+...
                    B(i,j)*T(i,j+1)+...
                    C(i,j)*T_next(i-1,j)+...
                    D(i,j)*T_next(i,j-1));
            end
        end
    end
    T_next(2:m-1,1)=T_next(2:m-1,2); % Left matrix column, Neumann BC
    T_next(2:m-1,n)=T_next(2:m-1,n-1); % Right matrix column, Neumann BC
    %r_iter = norm(T_next-T,2);
    r_iter = max(max(abs(T_next-T)));
    T = T_next;
    iter = iter+1;
end
fprintf('iterations: %4.1f \n',iter);

%% Plot Figures
figure
contourf(T)
title(['Gauss-Seidel Method, h =',num2str(h),', iterations: ',num2str(iter)])
xlabel('x points')
ylabel('y points')
colormap hot
colorbar('location','southoutside')