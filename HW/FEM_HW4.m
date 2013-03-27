%% FEM_HW4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving homework Problem No. 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% system k*u = f
k_activeDofs = [6 -3;-3 6];
f = [3/(2*pi^2);-3/(2*pi^2)];

% Solve for u
u_Dofs = k_activeDofs\f;

% x domain
x = 0:0.01:1;

% Compute elements displacements u(x)
u1 = x/(2*pi^2); u2 = (1-2*x)/(2*pi^2); u3 = (x-1)/(2*pi^2);
u = u1.*(x<1/3) + u2.*(x>1/3 & x<2/3) + u3.*(x>2/3);

% Exact solution
u_exact = 1/(pi^2)*(cos(pi*x)+2*x-1);

% Compute elements velocities du(x)
du1 = 1/(2*pi^2); du2 = -1/(pi^2); du3 = 1/(2*pi^2);
du = du1.*(x<1/3) + du2.*(x>1/3 & x<2/3) + du3.*(x>2/3);

% Exact Derivate of solution
du_exact = 1/(pi^2)*(2-pi*sin(pi*x));

% plot FEM solution and exact solution
figure(1); title('Solve ODE using 3 linear elements');
subplot(1,2,1); plot(x,u,'-dr',x,u_exact,'-.b'); ylabel('u'); xlabel('x'); legend('FEM','Exact')
subplot(1,2,2); plot(x,du,'-dr',x,du_exact,'-.b'); ylabel('du'); xlabel('x'); legend('FEM','Exact')


