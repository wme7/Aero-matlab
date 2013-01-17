function leg = legendre(x, p)

if(p==0)
    leg = ones(size(x));
end
if(p==1)
    leg = x;
end

if(p>=2) 
    legm1 = ones(size(x));
    leg   = x;
    
    for m=1:p-1
        legp1 = ((2*m+1)*x.*leg - m*legm1)/(m+1); 
        legm1 = leg;
        leg   = legp1;
    end
end

      