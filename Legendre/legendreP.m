function legP = legendreP(x,k)
% compute the value of the legendre polinomial of degree 'k' por column
% values 'x'

switch k
    case {0} % when k = 0
        legP = ones(size(x));
    case {1} % when k = 1
        legP = x;
    otherwise % k >= 2 
        legP0 = ones(size(x));
        legP1 = x;
        for m = (2):(k)
            legPX = ((2*m+1).*x.*legP1 - m.*legP0)/(m+1);
            legP0 = legP1;
            legP1 = legPX;
        end
        legP = legPX;
end