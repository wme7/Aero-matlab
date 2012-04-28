%% Quadrature
% Using quad to integrate f(x)

clear
clc

%% Evaluating quad
n=0;
for i=[2,5.5,9]
    n=n+1;
    F=@(x)((1./x)-1).^(i-1).*exp(1-(1./x))./(x.^2);
    [q(n,1) p(n,1)]=quad(F,0,1,1.0e-8);
end
[q p]

