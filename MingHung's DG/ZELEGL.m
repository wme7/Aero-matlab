function [ET,VN] = ZELEGL (N) 

% computes the nodes relative to the Legendre Gauss - Lobatto formula
% N  = order of the formula 
% ET = vector of the nodes , ET(I) , I = 1:N+1
% VN = values of the Legendre polynomial at the nodes , VN(I) , I = 1:N+1

if N == 0
    return ; 
end

N2 = fix ((N-1)/2) ;
SN = 2*N-4*N2-3 ;
ET(1) = -1 ;
ET(N+1) = 1 ;
VN(1) = SN ;
VN(N+1) = 1 ;

if N == 1
    return ;
end

ET(N2+2) = 0 ;
X = 0 ;
[Y,DY,D2Y] = VALEPO (N,X) ;
VN(N2+2) = Y ;
    
if N == 2
    return ;
end

C = pi/N ;
for i = 1:N2
    ETX = cos(C*i) ;
    for it = 1:8 
        [Y,DY,D2Y] = VALEPO (N,ETX) ;
        ETX = ETX-DY/D2Y ;
    end
    ET(i+1) = -ETX ;
    ET(N-i+1) = ETX ;
    VN(i+1) = Y*SN ;
    VN(N-i+1) = Y ;
end

