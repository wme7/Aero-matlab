function [x,w] = GaussLaguerre(n, alpha)

% This function determines the abscisas (x) and weights (w) for the
% Gauss-Laguerre quadrature of order n>1, on the interval [0, +infinity].
    % Unlike the function 'GaussLaguerre', this function is valid for
    % n>=34. This is due to the fact that the companion matrix (of the n'th
    % degree Laguerre polynomial) is now constructed as a symmetrical
    % matrix, guaranteeing that all the eigenvalues (roots) will be real.
    
    
% © Geert Van Damme
% geert@vandamme-iliano.be
% February 21, 2010    



% Building the companion matrix CM
    % CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial
    % under consideration. Moreover, CM will be constructed in such a way
    % that it is symmetrical.
i   = 1:n;
a   = (2*i-1) + alpha;
b   = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
CM  = diag(a) + diag(b,1) + diag(b,-1);

% Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
[V L]   = eig(CM);
[x ind] = sort(diag(L));
V       = V(:,ind)';
w       = gamma(alpha+1) .* V(:,1).^2;