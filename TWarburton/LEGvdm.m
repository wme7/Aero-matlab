function V = legVDM(x,p)

V = zeros(length(x), p+1);

for i=0:p
    V(:,i+1) = LEGENDRE(x,i);
end
