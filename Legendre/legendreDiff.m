function Dhat = legendreDiff(k)
% compute the coefficient Differentiation Matrix, using legendre polynomial
% of order 'k' as the base for our interpolation scheeme.

Dhat = zeros(k+1); % square matrix 

for j = 0:k
    for m = (j+1):2:k
        Dhat(j+1,m+1) = (2*j+1);
    end
end
        
