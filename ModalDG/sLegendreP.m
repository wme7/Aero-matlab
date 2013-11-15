function P = sLegendreP(k,x)
% Evaluate Legendre Polynomials of degree k, at a point x.

switch k
    case{0}
        P = 1;
    
    case{1}
        P = x;
        
    case{2}
        P = x.^2 - 1/12;
        
    case{3}
        P = x.^3 - 3/20*x;
        
    case{4}
        P = x.^4 - 3/14*x.^2 + 3/560;
        
    case{5}
        P = x.^5 - 5/18*x.^3 + 5/336*x;
        
    case{6}
        P = x.^6 - 15/44*x.^4 + 5/176*x.^2 - 5/14784;
        
    case{7}
        P = x.^7 -21/52*x.^5 + 105/2288*x.^3 - 35/27456*x;
        
    otherwise
        error('Out of range, max(k) = 10');
end
