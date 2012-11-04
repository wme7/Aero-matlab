function [PN,PD] = LPN (N,X)

PN(1) = 1 ;
PN(2) = X ;
PD(1) = 0 ;
PD(2) = 1 ;

P0 = 1 ;
P1 = X ;

for k = 2:N
    PF = (2*k-1)/k*X*P1-(k-1)/k*P0 ;
    PN(k+1) = PF ;
    if ( abs(X) == 1 )
        PD(k+1) = 0.5*X^(k+1)*k*(k+1) ;
    else
        PD(k+1) = k*(P1-X*PF)/(1-X*X) ;
    end
    P0 = P1 ;
    P1 = PF ;
end

