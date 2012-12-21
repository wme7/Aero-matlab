function [x,w] = sGaussLobatto(n)
%% Compute abscissas and weights Gauss-Lobatto Quadrature with n points
x = zeros(n,1); w = zeros(n,1); 
switch n 
    case{1} % 4th order Gauss-Lobatto quadrature
        x(1) = 0;
        w(1) = 1;
    case{2} % 4th order Gauss-Lobatto quadrature
        x(1) = -1/2;
        x(2) = 1/2;
        w(1) = 1/2;
        w(2) = 1/2;
    case{3} % 4th order Gauss-Lobatto quadrature
        % Abscissas
        x(1)=-1/2;
        x(2)= 0;
        x(3)= 1/2;
        % coefficients
        w(1)= 1/6;
        w(2)= 2/3;
        w(3)= 1/6;
    case{4} % 6th order Gauss-Lobatto quadrature
        % Abscissas
        x(1)=-0.5;
        x(2)=-sqrt(5)/10.0;
        x(3)= sqrt(5)/10.0;
        x(4)= 0.5;
        % coefficients
        w(1)= 1/12;
        w(2)= 5/12;
        w(3)= 5/12;
        w(4)= 1/12;
    case{5} % 8th order Gauss-Lobatto quadrature
        % Abscissas
        x(1)=-0.5;
        x(2)=-sqrt(21)/14;
        x(3)= 0;
        x(4)= sqrt(21)/14;
        x(5)= 0.5;
        % coefficients
        w(1)= 1/20;
        w(2)= 49/180;
        w(3)= 64/180;
        w(4)= 49/180;
        w(5)= 1/20;
    case{6} % 10th order Gauss-Lobatto quadrature
        % Abscissas
        x(1)=-1/2;
        x(2)=-sqrt(147+42*sqrt(7))/42;
        x(3)=-sqrt(147-42*sqrt(7))/42;
        x(4)= sqrt(147-42*sqrt(7))/42;
        x(5)= sqrt(147+42*sqrt(7))/42;
        x(6)= 1/2;
        % coefficients
        w(1)= 1/30;
        w(2)= (7 +5*sqrt(7))*sqrt(7)/(7+sqrt(7))/20;
        w(3)= (-7+5*sqrt(7))*sqrt(7)*(7+sqrt(7))/840;
        w(4)= (-7+5*sqrt(7))*sqrt(7)*(7+sqrt(7))/840;
        w(5)= (7 +5*sqrt(7))*sqrt(7)/(7+sqrt(7))/20;
        w(6)= 1/30;
    otherwise
        error('out of range, max(n) = 6')
end