%% FEM class HW1
% In this subroutine 1D stressed beam and axisymmetric heated plate.
%
% Coded by Manuel Diaz, f99543083, IAM, 2013.03.05.
clear all; close all; clc;

%% Exercise No.1
x = linspace(1,3,40); % Domain
A = 0.01;    % Cross section Area
E = 10^5;    % Young Modulus

% Exact Solution
c_1 = 1/1000 + 9/(A*E) - 3*(A+10)/(10*A*E);
c_2 = (A+10)/(10*A*E);
u_exact = c_2*x + c_1 - x.^3/(3*A*E);
s_exact = E*c_2 - x.^2/A; 

% Linear Approximation
alpha1_1 = (3*A-100)/(30*A*E);
u1 = 1/1000 + alpha1_1*(x-3);
s1 = E*(alpha1_1*ones(size(x)));

% Parabolic Approximation
alpha2_1 = -(220-3*A)/(30*A*E);
alpha2_2 = -2/(A*E);
u2 = 1/1000 + alpha2_1*(x-3) + alpha2_2*(x-3).^2;
s2 = E*(alpha2_1 + 2*alpha2_2*(x-3));

% Draw Figures
figure(1)
subplot(1,2,1); 
hold on; title('Displacement, u(x)')
plot(x,u_exact,'-k'); plot(x,u1,'ok'); plot(x,u2,'sk'); hold off;
xlabel('x'), ylabel('displacement'); legend('exact','linear','quadratic',3)
subplot(1,2,2); 
hold on; title('Stress, \sigma(x)')
plot(x,s_exact,'-k'); plot(x,s1,'ok'); plot(x,s2,'sk'); hold off;
xlabel('x'), ylabel('stress'); legend('exact','linear','quadratic',3)

%% Exercise No.2
R = 2/1000;  % Radious of plate [m]
r = linspace(0,R,40); % Domain
s = 3.18E9;  % Source term [W/m^3]
k = 15;      % Conduction constant [W/(m.K)]
C_1 = 0;
C_2 = s*R^2/(2*(k+1));
T_exact = C_1*r.^(1-1/k)+C_2-s*r.^2/(2*(k+1));

% Linear Approximation
Alpha1_1 = s/(750*k);
Alpha1_0 = - Alpha1_1*(2/1000);
T1 = Alpha1_0 + Alpha1_1*r;

% Quadratic Approximation
 Alpha2_1 = 0.006*s/k;
Alpha2_2 = - 1.75*s/k;
Alpha2_0 = - Alpha2_1*(2/1000) - Alpha2_2*(2/1000)^2;
T2 = Alpha2_0 + Alpha2_1*r + Alpha2_2*r.^2;

% Draw Figure
figure(2)
hold on; title('Plate Temperature, T(r)')
plot(r,T_exact,'-k'); plot(r,T1,'ok'); plot(r,T2,'sk');
xlabel('r'),ylabel('temperature ºK');legend('exact','linear','quadratic',4)