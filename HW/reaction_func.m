function dc=reaction(t,c)
%% Convencion
% c(1) = Cs; c(2) = Cp, c(3) = Ces, c(4) = Cp
%
%% Constants
k1 = 2.0e3; k2 = 5.5e-5; k3 = 10.0;
a=[-k1 -k1 k1 0]; b=[k2 (k2+k3) -(k2+k3) k3];
%% Function
dc = zeros(4,1);
dc(1) = a(1)*c(1)*c(2) + b(1)*c(3);
dc(2) = a(2)*c(1)*c(2) + b(2)*c(3);
dc(3) = a(3)*c(1)*c(2) + b(3)*c(3);
dc(4) = a(4)*c(1)*c(2) + b(4)*c(3);
