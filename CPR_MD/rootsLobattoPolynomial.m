function xi = rootsLobattoPolynomial(k_deg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find roots of Lobatto Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input : k_deg: Lobatto Polynomial     
% Output : xi: roots of Lobatto polynomial, Lobatto Solution Points
%    

switch k_deg
    case 2
    lobattoP2 = 3/2*[1 0 -1]; xi = sort(roots(lobattoP2));
    case 3
    lobattoP3 = 5/2*[1 0 -1 0]; xi = sort(roots(lobattoP3));
    case 4
    lobattoP4 = 7/8*[5 0 -6 0 1]; xi = sort(roots(lobattoP4));
    case 5
    lobattoP5 = 9/8*[7 0 -10 0 3]; xi = sort(roots(lobattoP5));
    otherwise 
        error('Lobatto Polynomial case not available')
end
%% Plot Lobatto Polynomial (test*)
%Interval = -1:0.1:1; plot(Interval,polyval(xi,Interval)); grid on;

    
    
    