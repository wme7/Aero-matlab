function [WT] = WELEGL (N,ET,VN)

% computes the weights relative to the Legendre Gauss - Lobatto formula
% N  = order of the formula 
% ET = Jacobi Gauss - Lobatto nodes , ET(I) , I = 1:N+1
% VN = values of the Legendre polynomial at the nodes , VN(I) , I = 1:N+1
% WT = vector of the weights , WT(I) , I = 1:N+1

if N == 0
    return ;
end

N2 = fix ((N-1)/2) ; 
DN = N ;
C = 2/(DN*(DN+1)) ;
for i = 0:N2
    X = ET(i+1);
    Y = VN(i+1);
    WTX = C/(Y*Y) ;
    WT(i+1) = WTX ;
    WT(N-i+1) = WTX ;
end

if (N-1) == (2*N2)
    return ;
end

X = 0 ;
Y = VN(N2+2) ;
WT(N2+2) = C/(Y*Y) ;