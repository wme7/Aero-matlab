function P = LegMat(k,x)
% Construct Array P
%**************************************************************************
% P is a matrix array in which each element is a scaled Legendre
% Polynomials of Degree (l=0:k) evaluated at a column array of points (x)
%
% Coded by Manuel Diaz 2012.12.05
%**************************************************************************
%% Allocate Solution Variable
n = length(x); P = zeros(n,n);

%% Evaluate sLegendre fucntion
for l = 0:k         % all Polynomial Degrees up to k
    j = l + 1;      % Dummy index
    for i = 1:n     % all Element Points in x
        P(j,i) = legendreP(x(i),l);
    end
end
P = P';             % Transpose our result