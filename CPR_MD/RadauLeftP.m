function RLP = RadauLeftP(kDeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load from table the coefs for Radau Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : kDeg: Polynomial Degree requested     
% Output : RLP: Left Radau polynomial
%

RLP = (1/2)*( LegendreP(kDeg) + LegendreP(kDeg-1) );

%% Plot Polynomial (test*)
%ezplot(RLP,[-1,1]); grid on;  