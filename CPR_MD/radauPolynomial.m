function radauP = radauPolynomial(k_deg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load from table the coefs for Radau Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : k_deg: Radau Polynomial     
% Output : radauP: Radau polynomial coeficients
%    

switch k_deg
    case(1)
    radauP = 1/2*[1 1];
    case(2)
    radauP = 1/4*[3 2 -1];
    case(3)
    radauP = 1/4*[5 3 -3 -1];
    case(4)
    radauP = 1/16*[35 20 -30 -12 3];
    otherwise
        error('Radau Polynomial coefcients not available')
end    
%% Plot Radau Polynomial (test*)
%Interval = -1:0.1:1; plot(Interval,polyval(radauP,Interval)); grid on;   