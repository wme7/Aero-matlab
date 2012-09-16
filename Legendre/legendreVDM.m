function V = legendreVDM(x,k)
% Compute the vandermonde matrix of order k for the legendre polinomials of
% a column vector of modal points 'x'
%
% input:  row vector x \subset (-1,1) of points k+1, where k is 
%         the maximal order of the Legendre polynomials
% output: a matrix P = P(length(x),k+1) containing the values of the
%         Legendre polynomials at the points x
%

V = zeros(length(x),k+1);
for m = 0:k
    V(:,m+1) = legendreP(x,m);
end