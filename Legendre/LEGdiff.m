function  Dhat = LEGdiff(p)

Dhat = zeros(p+1);

for n=0:p
    for m=n+1:2:p
          Dhat(n+1,m+1) =  (2*n+1);
    end
end
      