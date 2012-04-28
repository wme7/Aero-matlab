%% 2D Poission Equation
% This program solves the following PDE:
%
% $$\frac{\partial^2{T}}{\partial{x^2}}+\frac{\partial^2{T}}{\partial{y^2}}=0$$

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

% initaliza Matriz T
T=zeros(n+1,m+1);

%% Boundary Contidions
% Our domain requires 4 boudary conditions in total, can be only nuemann, 
% or only dirichlet or a combination of both.

%% Nuemann Boundary Conditions (Prescribed values at the borders)
% Bottom (last row of the matrix)
T(1,:)=0;
% Top (first row of the matrix)
%T(1,:)=1;
for i=0:n
    %using function of x to prescribe the values.
    T(m+1,i+1)=1-cos(pi/(n)*i); 
end
% Left
%T(:,1)=0;
% Right
%T(:,n+1)=2;

%T

% Vector of Prescribed Bondary Values
% Bottom
fb=(T(1,2:n)).' ; %transpose this array
% Top
ft=T(m+1,2:n).' ; %transpose this array
% Left
fl=T(2:n,1);
% Right
fr=T(2:n,n+1);

% Computing vector of known values
fs=zeros(n-1,1);
fs(1)=fl(1);
fs(n-1)=fr(1);
w=fb+fs;
for i=2:n-2
    fs(1)=fl(i);
    fs(n-1)=fr(i);
    w=cat(1,w,fs);
end
fs=zeros(n-1,1);
fs(1)=fl(n-1);
fs(n-1)=fr(n-1);
y=ft+fs;
w=cat(1,w,y);

%% Dirichlect Boundary Conditions (flux at the borders of the domain)
% This conditions requires de mofification of the diagonal values of our
% solution matrix:

% Bottom (last row of the matrix)
%sb=ones(n-1,1); % to insulate this boundary
sb=zeros(n-1,1); %to suppres this condition
% Top (first row of the matrix)
%st=ones(n-1,1); % to insulate this boundary
st=zeros(n-1,1); %to suppres this condition
% Left
sl=zeros(m-1,1); 
sl(1)=1; % to insulate this boundary
% Right
sr=zeros(m-1,1);
sr(m-1)=1; % to insulate this boundary

% Computing the diagonal matrix to substract
s=sl+sr; %sum condition in the right and left
stt=st+s; %sum condition of the top and sides
sbb=sb+s; %sum condition of the bottom and sides
u=cat(1,s,s);
for i=4:n-2
    u=cat(1,u,s);
end
z=cat(1,stt,u);
z=cat(1,z,sbb);
Z=spdiags(z,0,length(z),length(z)); %transforms the vector in a diagonal mat.
%spy(Z)

%% System Solution
% Our solution has the shape: 
%
% $$ A_{ij}u_j=f_i  $$

% Formulating A:
I=eye(m-1);
e=ones(m-1,1);
K=spdiags([e,-4*e,e],[-1,0,1],m-1,m-1);
S=spdiags([e,e],[-1,1],m-1,m-1);
A=(kron(I,K)+kron(S,I)); 
A=(A+Z); %(for flux conditions)
u=inv(A)*(-w);
spy(A)

%% Re-arranging the data in vector u

for j=1:n-1
    for i=1:n-1
        t(j,i)=u(i+(n-1)*(j-1));
    end
end
for j=1:n-1
    for i=1:n-1
        T(i+1,j+1)=t(i,j);
    end
end

% Make up for the plot
if sum(sb)>=1
    T(1,:)=T(2,:);
end
if sum(st)>=1
    T(n+1,:)=T(n,:);
end
if sum(sl)>=1
    T(:,1)=T(:,2);
end
if sum(sr)>=1
    T(:,n+1)=T(:,n);
end
%T
%surf(T)
contourf(T)
colormap hot
colorbar('location','southoutside')

%% Computing the Error of our aproximation
% % The exact solution of our problem is given by:
% %
% % $$ T(x,y)=y-\frac{sinh(\pi y))}{sinh(\pi)}cos(\pi x) $$
% for j=1:m+1
%     for i=1:n+1
%         T2(j,i)=j/(m+1)-sinh(pi/(m+1)*j)/sinh(pi)*cos(pi/(n+1)*i);
%     end
% end
% Error=T2-T;
% %T2
% h
% NormError=norm(Error)
% %surf(T2)
% contourf(T2)
% colormap hot
% colorbar('location','southoutside')


