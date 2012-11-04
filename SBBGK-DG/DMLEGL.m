function [DMA] = DMLEGL (N,NM,ET,VN)

% computes the entries of the derivative matrix relative to the 
% Legendre Gauss - Lobatto nodes
% N   = parameter relative to the dimension of the matrix
% NM  = order of the matrix as declared in the main dimension statement
% ET  = vector of the nodes , ET(I) , I = 1:N+1
% VN  = values of the Legendre polynomial at the nodes , VN(I) , I = 1:N+1
% DMA = derivative matrix , DMA(I,J) , I = 1:N+1 , J = 1:N+1

DMA(1,1) = 0 ;

if N == 0
    return ;
end

for i = 0:N
    VI = VN(i+1) ;
    EI = ET(i+1) ;
    for j = 0:N
        if i ~= j
            VJ = VN(j+1) ;
            EJ = ET(j+1) ;
            DMA(i+1,j+1) = VI/(VJ*(EI-EJ)) ;
        else
            DMA(i+1,i+1) = 0 ;
        end
    end    
end

DN = N ;
C = 0.25*DN*(DN+1) ;
DMA(1,1) = -C ;
DMA(N+1,N+1) = C ;

