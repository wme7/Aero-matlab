function lobP = LobattoP(kDeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lobatto Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : kDeg: Polynomial Degree requested     
% Output : lobP: Symbolic Lobatto polynomial
%

lobP = LegendreP(kDeg) - LegendreP(kDeg-2);

%% Plot Polynomial (test*)
ezplot(lobP,[-1,1]); grid on; 