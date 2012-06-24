function udot=LorenzAttractorfunc(t,u)
% Lorenz Attractor system of equations.
%
% Convention: u(1) = x, u(2) = y, u(3) = z 
%
sigma = u(4); b = u(5); r = u(6);

udot=zeros(6,1); % !because our imput is a 6x1 row
udot(1)= sigma*(u(2)-u(1));
udot(2)= u(1)*(r-u(3))-u(2);
udot(3)= u(1)*u(2)-b*u(3);