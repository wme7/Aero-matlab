
% purpose: solve du/dt = mu*u
%                u(0)  = uo
%          for t in [0,T] with intervals of dt   

% input: 
% T  = final time
% dt = time step 
% uo = initial value of solution
% mu = multiplier for ODE rhs
% output:
% u  = vector of values (u_0, u_1,... , u_nsteps)
% dt = actual time step used

function [u,t] = EulerForwardODE(T, dt, uo, mu)

 nsteps = ceil(T/dt);
 dt     = T/nsteps;
 u      = zeros(nsteps+1,1);
 t      = zeros(nsteps+1,1);
 
 u(1) = uo;
 t(1) = 0;
 for n=1:nsteps
     f = mu*u(n);
     u(n+1) = u(n) + dt*f;
     t(n+1) = t(n) + dt;
 end
 
 plot(t,abs(u), 'k-*');
 hold on;
 plot(t,real(u), 'r-*'); 
 plot(t,imag(u), 'b-*'); 
 plot(t,abs(exp(t*mu)), 'g-*'); % exact solution 
 hold off;
 legend('|u_n|', 'Re(u_n)', 'Im(u_n)', '|u|');

 xlabel('t', 'FontSize', 18);
 tt = sprintf('Euler-Forward: dt=%g, mu=%g+%g*i', ...
	     dt, real(mu), imag(mu));
 title(tt)
