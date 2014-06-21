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
 J(i,j) = xn(i,j)*ye(i,j)-yn(i,j)*xe(i,j); % Jacobian
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
view(150,30)
title('Jacobian ','FontSize',18);