function chebP = chebyshevP(x,k)
% compute the value of the chebyshev polinomial of degree 'k' for every
% column value 'x'

%number of points
p = k+1;

switch k
    case {0} % when k = 0
        chebP = ones(size(x));
    case {1} % when k = 1
        chebP = x;
    otherwise % k >= 2 
        chebP0 = ones(size(x));
        chebP1 = x;
        for m = 1:k-1
            chebPX = (2*x.*chebP1 - chebP0);
            chebP0 = chebP1;
            chebP1 = chebPX;
        end
        chebP = chebPX;
end