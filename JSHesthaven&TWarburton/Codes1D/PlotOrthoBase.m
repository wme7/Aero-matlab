% plot Jacobi Polynomial in 2D

r = (-1:0.25:1)'; s = (-1:0.25:1)'; %[r,s] = meshgrid(s,r);

a = 2*(1+r)./(1-s)-1; b = s;
a(end)=-1;

i=1; j=1;

h1 = JacobiP(a,0,0,i); h2 = JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;

scatter3(a,b,P);