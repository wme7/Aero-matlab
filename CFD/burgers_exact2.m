function burgers_exact2

%SIDENOTE
%If the characteristics of the exact solution do not intersect, then the
%classical solution of burgers equation is given by:

% Domain
x = 0:0.02:1;
% IC (at t=0)
u0 = 0.2 + sin(2*pi*x);
% u function
u = @(x,t) 0.2 + sin(2*pi*x);

for t = 0:0.01:0.1
    u_next = u(x - t*u(x,t),t);
    plot(x,u_next,'o');
    pause(0.1)
    drawnow
end

% This solution is valid for any u(x,t) function.
% However if the characteristics colide, then the classical solution is not
% valid anymore. It can only be solved computationally.