function xdot=dpendulum(t,x)
%% Double Pendulum ODE
%
% ivp=[gamma0; dtgamma0; alpha0; dtalpha0; beta0; lambda; omega; psi; eta];

beta0=x(5); lambda=x(6); omega=x(7); psi=x(8); eta=x(9);
xdot=zeros(9,1); % a column vector

C = cos(beta0)*cos(x(3)-x(1))-sin(beta0)*sin(x(3)-x(1));
S = sin(beta0)*cos(x(3)-x(1))+cos(beta0)*sin(x(3)-x(1));

xdot(1) = x(2);
xdot(2) = ((omega^2)*sin(x(1))-psi*(S*(x(4)^2+C*x(2)^2*eta)+...
    C*(lambda^2)*sin(x(3))))/(-1+(C^2)*eta*psi);
xdot(3) = x(4);
xdot(4) = (S*eta*(x(2)^2+C*x(4)^2*psi)-C*eta*(omega^2)*sin(x(1))+...
    (lambda^2)*sin(x(3)))/(-1+(C^2)*eta*psi);