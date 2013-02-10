nx = length(x);
dx = max(abs(x(1:end-1)-x(2:end)));
cfl = 0.6;
a = -1;
dt = dx*cfl/abs(a);
un = zeros(size(u));

eps = 1E-25;
c1 = - 1/6;
c2 = 2/6;
c3 = 5/6;
c4 = - 7/6;
c5 = 11/6;

k = 1; %:nv
    
vxp = max(a,0);
vxm = min(a,0);

for j = 4:nx-3
    
sl0m = vxm*(13./12.)*(u(k,j-2)-2*u(k,j-1)+u(k,j))^2 + vxm*(1./4.)*(u(k,j-2)-4.*u(k,j-1)+3.*u(k,j))^2;
sl1m = vxm*(13./12.)*(u(k,j-1)-2*u(k,j)+u(k,j+1))^2 + vxm*(1./4.)*(u(k,j-1)-u(k,j+1))^2;
sl2m = vxm*(13./12.)*(u(k,j)-2*u(k,j+1)+u(k,j+2))^2 + vxm*(1./4.)*(3.*u(k,j)-4.*u(k,j+1)+u(k,j+2))^2;
sl0p = vxp*(13./12.)*(u(k,j-3)-2*u(k,j-2)+u(k,j-1))^2 + vxp*(1./4.)*(u(k,j-3)-4.*u(k,j-2)+3.*u(k,j-1))^2;
sl1p = vxp*(13./12.)*(u(k,j-2)-2*u(k,j-1)+u(k,j))^2 + vxp*(1./4.)*(u(k,j-2)-u(k,j))^2;
sl2p = vxp*(13./12.)*(u(k,j-1)-2*u(k,j)+u(k,j+1))^2 + vxp*(1./4.)*(3.*u(k,j-1)-4.*u(k,j)+u(k,j+1))^2;

sr0m = vxm*(13./12.)*(u(k,j-1)-2*u(k,j)+u(k,j+1))^2 + vxm*(1./4.)*(u(k,j-1)-4.*u(k,j)+3.*u(k,j+1))^2;
sr1m = vxm*(13./12.)*(u(k,j)-2*u(k,j+1)+u(k,j+2))^2 + vxm*(1./4.)*(u(k,j)-u(k,j+2))^2;
sr2m = vxm*(13./12.)*(u(k,j+1)-2*u(k,j+2)+u(k,j+3))^2 + vxm*(1./4.)*(3.*u(k,j+1)-4.*u(k,j+2)+u(k,j+3))^2;
sr0p = vxp*(13./12.)*(u(k,j-2)-2*u(k,j-1)+u(k,j))^2 + vxp*(1./4.)*(u(k,j-2)-4.*u(k,j-1)+3.*u(k,j))^2;
sr1p = vxp*(13./12.)*(u(k,j-1)-2*u(k,j)+u(k,j+1))^2 + vxp*(1./4.)*(u(k,j-1)-u(k,j+1))^2;
sr2p = vxp*(13./12.)*(u(k,j)-2*u(k,j+1)+u(k,j+2))^2 + vxp*(1./4.)*(3.*u(k,j)-4.*u(k,j+1)+u(k,j+2))^2;

al0m = 1 / (10 * (eps + sl0m))^2;
al1m = 6 / (10 * (eps + sl1m))^2;
al2m = 3 / (10 * (eps + sl2m))^2;
al0p = 3 / (10 * (eps + sl0p))^2;
al1p = 6 / (10 * (eps + sl1p))^2;
al2p = 1 / (10 * (eps + sl2p))^2;

ar0m = 1 / (10 * (eps + sr0m))^2;
ar1m = 6 / (10 * (eps + sr1m))^2;
ar2m = 3 / (10 * (eps + sr2m))^2;
ar0p = 3 / (10 * (eps + sr0p))^2;
ar1p = 6 / (10 * (eps + sr1p))^2;
ar2p = 1 / (10 * (eps + sr2p))^2;

wl0m = al0m / (al0m+al1m+al2m);
wl1m = al1m / (al0m+al1m+al2m);
wl2m = al2m / (al0m+al1m+al2m);
wl0p = al0p / (al0p+al1p+al2p);
wl1p = al1p / (al0p+al1p+al2p);
wl2p = al2p / (al0p+al1p+al2p);

w0(j) = wl0m;
w1(j) = wl1m;
w2(j) = wl2m;

wr0m = ar0m / (ar0m+ar1m+ar2m);
wr1m = ar1m / (ar0m+ar1m+ar2m);
wr2m = ar2m / (ar0m+ar1m+ar2m);
wr0p = ar0p / (ar0p+ar1p+ar2p);
wr1p = ar1p / (ar0p+ar1p+ar2p);
wr2p = ar2p / (ar0p+ar1p+ar2p);

ulm(j) = vxm *(wl0m*(c1*u(k,j-2)+c3*u(k,j-1)+c2*u(k, j ))+...
               wl1m*(c2*u(k,j-1)+c3*u(k, j )+c1*u(k,j+1))+...
               wl2m*(c5*u(k, j )+c4*u(k,j+1)+c2*u(k,j+2)));
ulp(j) = vxp *(wl0p*(c2*u(k,j-3)+c4*u(k,j-2)+c5*u(k,j-1))+...
               wl1p*(c1*u(k,j-2)+c3*u(k,j-1)+c2*u(k, j ))+...
               wl2p*(c2*u(k,j-1)+c3*u(k, j )+c1*u(k,j+1)));
urm(j) = vxm *(wr0m*(c1*u(k,j-1)+c3*u(k, j )+c2*u(k,j+1))+...
               wr1m*(c2*u(k, j )+c3*u(k,j+1)+c1*u(k,j+2))+...
               wr2m*(c5*u(k,j+1)+c4*u(k,j+2)+c2*u(k,j+3)));
urp(j) = vxp *(wr0p*(c2*u(k,j-2)+c4*u(k,j-1)+c5*u(k, j ))+...
               wr1p*(c1*u(k,j-1)+c3*u(k, j )+c2*u(k,j+1))+...
               wr2p*(c2*u(k, j )+c3*u(k,j+1)+c1*u(k,j+2)));

ul(j) = ulm(j) + ulp(j);
ur(j) = urm(j) + urp(j);

un(j) = u(k,j) - dt/dx * (ur(j) - ul(j));

end