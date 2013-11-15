function dP = dsLegendreP(k,x)
% Evaluate the Derivative of our scaled Legendre Polynomials of order k, at a point x.

switch k
    case{0}
        dP = zeros(length(x),1);
    
    case{1}
        dP = ones(length(x),1);
        
    case{2}
        dP = 2*x;
        
    case{3}
        dP = 3*x.^2-3/20;
        
    case{4}
        dP = 4*x.^3 - 3/7*x;
        
    case{5}
        dP = 5*x.^4 - 5/3*x.^2 + 5/336;
        
    case{6}
        dP = 6*x.^5 - 15/11*x.^3 + 5/176*2*x;
        
    case{7}
        dP = 7*x.^6 -21/52*5*x.^4 + 105/2288*3*x.^2 - 35/27456;
        
    otherwise
        error('Out of range, max(k) = 10');
end
