function RRP = RadauRightP(kDeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load from table the coefs for Radau Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : kDeg: Polynomial Degree requested
% Output : RRP: Right Radau polynomial
%

RRP = (-1)^(kDeg)/2*( LegendreP(kDeg) - LegendreP(kDeg-1) );

%% Plot Polynomial (test*)
ezplot(RRP,[-1,1]); grid on;  