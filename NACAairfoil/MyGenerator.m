close all;
clear all;
 
%-------------------------------------------------------------------------% 
% Define all the parameters 
%-------------------------------------------------------------------------% 
 
Na = 50;
M = 30;
A = 10; %Clustering points near leading edge, inf = no clustering
B = 2; %Clustering points near the surface, 1 = no clustering
t = .12;
c = 1;
L = 1;
N = floor(2*c/(2*c+L)*Na);
x1 = c*ones(1,N+1);
y1 = zeros(1,N+1);
x2 = x1;
y2 = y1;
y2(1) = c;
s1 = y1;
s2 = s1;
N2 = 250;

%-------------------------------------------------------------------------% 
% Solve Algebraic Method
%-------------------------------------------------------------------------%
for i = 2:N2+1
    x1(i) = c*(1-(i-1)/N2);
    y1(i) = max([t/.2*(.2969*x1(i)^.5-.126*x1(i)-.3516*x1(i)^2+.2843*x1(i)^3-.1015*x1(i)^4) 0]);
    s1(i) = s1(i-1)+sqrt((x1(i)-x1(i-1))^2+(y1(i)-y1(i-1))^2);
    x2(i) = c*(1-2*(i-1)/N2);
    if x2(i) >= 0;
        y2(i) = c;
    else
        y2(i) = sqrt(c^2-(x2(i))^2);
    end
    s2(i) = s2(i-1)+sqrt((x2(i)-x2(i-1))^2+(y2(i)-y2(i-1))^2);
end
S1 = s1(end)*(1-exp(-(0:N)/A))/(1-exp(-(N)/A));
S1(end)/s1(end)
S2 = (0:N)*s2(end)/N;
for i = 1:N+1
    X1(i) = interp1(s1,x1,S1(i));
    Y1(i) = interp1(s1,y1,S1(i));
    X2(i) = interp1(s2,x2,S2(i));
    Y2(i) = interp1(s2,y2,S2(i));
end
r = exp(((1:(M+1))-M+1)/((M+1)/B));
r = r-r(1);
r = r/max(r);
for i = 1:N+1
    for j = 1:M+1
        x(i,j) = X1(i) + r(j)*(X2(i)-X1(i));
        y(i,j) = Y1(i) + r(j)*(Y2(i)-Y1(i));
    end
end
N = Na-N;
x = [((L+c):-L/N:(c+L/N))'*ones(1,M+1); x];
y = [(c*r'*ones(1,N))'; y];
x = [x; x(end-1:-1:1,:)];
y = [y ;-y(end-1:-1:1,:)];
[n,m] = size(x);

%-------------------------------------------------------------------------% 
% Plot Algebraic Method
%-------------------------------------------------------------------------%
mesh(x,y,zeros(n,m))
view(0,90)
colormap([0 0 0])
axis('equal','tight');
xlabel('Length [unitless]','FontSize',14);
ylabel('Length [unitless]','FontSize',14);
title('Algebraic Method ','FontSize',18);

%-------------------------------------------------------------------------% 
% Predetermined Constants
%-------------------------------------------------------------------------%
de = 1/M;
de2 = de^2;
dn = 1/Na;
dn2 = dn^2;
dedn = 2/(M*Na);
P = 100; %number of iterations for the elliptic solver.

%-------------------------------------------------------------------------% 
% Solve Elliptic Method
%-------------------------------------------------------------------------%
 
for k = 1:P
    for i = 2:2*Na-1
        for j = 2:M-1
            xn = (x(i+1,j)-x(i-1,j))/(2*dn);
            yn = (y(i+1,j)-y(i-1,j))/(2*dn);
            xe = (x(i,j+1)-x(i,j-1))/(2*de);
            ye = (y(i,j+1)-y(i,j-1))/(2*de);
            a(i,j) = xn^2+yn^2;
            b(i,j) = xe*xn+ye*yn;
            c(i,j) = xe^2+ye^2;
        end
    end
    for i = 2:2*Na-1
        for j = 2:M-1
            x(i,j) = (a(i,j)/de2*(x(i+1,j)+x(i-1,j))+c(i,j)/dn2*(x(i,j+1)+x(i,j-1))...
                -b(i,j)/dedn*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1)))...
                /(2*(a(i,j)/de2+c(i,j)/dn2));
            y(i,j) = (a(i,j)/de2*(y(i+1,j)+y(i-1,j))+c(i,j)/dn2*(y(i,j+1)+y(i,j-1))...
                -b(i,j)/dedn*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1)))...
                /(2*(a(i,j)/de2+c(i,j)/dn2));
        end
    end
end
xn = ones(Na,M)/0;
yn = xn;
xe = xn;
ye = xn;
for i = 2:Na-1
    for j = 2:M-1
        xn(i,j) = (x(i+1,j)-x(i-1,j))/(2*dn);
        yn(i,j) = (y(i+1,j)-y(i-1,j))/(2*dn);
        xe(i,j) = (x(i,j+1)-x(i,j-1))/(2*de);
        ye(i,j) = (y(i,j+1)-y(i,j-1))/(2*de);
        J(i,j) = yn(i,j)*xe(i,j)-xn(i,j)*ye(i,j); % Jacobian
    end
end
figure(2)
mesh(x,y,zeros(n,m))
xlabel('Length [unitless]','FontSize',14);
ylabel('Length [unitless]','FontSize',14);
title('Elliptic Method ','FontSize',18);
axis('equal','tight');

%-------------------------------------------------------------------------% 
% Plot Elliptic Method
%-------------------------------------------------------------------------%
 
view(0,90)
colormap([0 0 0])

%-------------------------------------------------------------------------% 
% Plot Jacobian
%-------------------------------------------------------------------------%
 
figure(3)
x=2:30; y=2:50; [X,Y]=meshgrid(x,y);
surf(X,Y,J);
view(-45,45)
title('Jacobian ','FontSize',18);
