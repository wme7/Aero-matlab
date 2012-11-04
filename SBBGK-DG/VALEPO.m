function [Y,DY,D2Y] = VALEPO (N,X)

% computes the value of the Legendre polynomial of degree N
% and its first and second derivatives at a given point
% N   = degree of the polynomial
% X   = point in which the computation is performed 
% Y   = value of the polynomial in X
% DY  = value of the first derivative in X
% D2Y = value of the second derivative in X

Y = 1;
DY = 0 ;
D2Y = 0 ;

if N == 0
    return ;
end

Y = X ;
DY = 1 ;
D2Y = 0 ;

if N == 1
    return ;
end

YP = 1 ;
DYP = 0 ;
D2YP = 0 ;

for i = 2:N
    C1 = i ;
    C2 = 2*C1-1 ;
    C4 = C1-1 ;
    YM = Y ;
    Y = (C2*X*Y-C4*YP)/C1 ;
    YP = YM ;
    DYM = DY ;
    DY = (C2*X*DY-C4*DYP+C2*YP)/C1 ;
    DYP = DYM ;
    D2YM = D2Y ;
    D2Y = (C2*X*D2Y-C4*D2YP+2*C2*DYP)/C1 ;
    D2YP = D2YM ;
end    
