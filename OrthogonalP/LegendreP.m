function legP = LegendreP(kDeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legendre Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : kDeg: Polynomial Degree requested     
% Output : legP: symbolic Legendre polynomial
%
x = sym('x');
switch kDeg
    case(0)
        legP = x^0; % = 1
    case(1)
        legP = x;
    case(2)
        legP = 1/2*(-1+3*x^2);
    case(3)
        legP = 1/2*(-3*x+5*x^3);
    case(4)
        legP = 1/8*(3-30*x^2+35*x^4);
    case(5)
        legP = 1/8*(63*x^5-70*x^3+15*x);
    case(6)
        legP = 1/16*(231*x^6-315*x^4+105*x^2-5);
    otherwise
        error('Legendre Polynomial not available');
end
%% Plot Polynomial (test*)
ezplot(legP,[-1,1]); grid on; 