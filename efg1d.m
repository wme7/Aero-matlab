% ONE DIMENSIONAL EFG PROGRAM
% SET UP NODAL COORDINATES ALONG BAR, DETERMINE NUMBER OF CELLS
x = [0.0:.1:1.0];
nnodes = length(x);         
ncells = nnodes-1;
nosDeDirichlet = 1;

% SET PARAMETERS FOR WEIGHT FUNCTION, MATERIAL PROPERITES
dmax = 2.0;
E=1.0; area=1.0;

% DETERMINE DMI FOR EACH NODE
dm = dmax*(x(2)-x(1))*ones(1,nnodes);

%SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
gg = zeros(1,ncells);
jac = (x(2)-x(1))/2;
weight = 2;
gg = -.05:.1:0.95; gg(1) = 0.0;

% INITIALIZE MATRICES
k = zeros(nnodes);
f = zeros(nnodes,1);
GG = zeros(nnodes,nosDeDirichlet);

% LOOP OVER GAUSS POINTS
for j = 1:length(gg)
   xg = gg(j);

% DETERMINE DISTANCE BETWEEN NODES AND GAUSS POINT
dif = xg*ones(1,nnodes)-x;
 
% SET UP WEIGHTS W AND DW FOR EACH NODE
clear w dw
for i=1:nnodes
drdx = sign(dif(i))/dm(i);
r = abs(dif(i))/dm(i);
if r<=0.5
  w(i) = (2/3) - 4*r*r + 4*r^3;
  dw(i) = (-8*r + 12*r^2)*drdx;
elseif r<=1.0
  w(i) = (4/3)-4*r+4*r*r -(4/3)*r^3;
  dw(i) = (-4 + 8*r-4*r^2)*drdx;
elseif r>1.0
w(i) = 0.0;
dw(i) = 0.0;
end
end

%SET UP SHAPE FUNCTIONS AND DERIVATIVES
won = ones(1,nnodes);
p = [won;x];
B = p.*[w;w];
pp = zeros(2);
A = zeros(2);
dA = zeros(2);
for i=1:nnodes
   pp = p(1:2,i)*p(1:2,i)';
   A = A+w(1,i)*pp;
   dA = dA+dw(1,i)*pp;
end
Ainv = inv(A);
pg = [1 xg];
phi(j,:) = pg*Ainv*B;
db = p.*[dw;dw];
da = -Ainv*(dA*Ainv);
dphi(j,:) = [0 1]*Ainv*B+pg*(da*B+Ainv*db);

%ASSEMBLE DISCRETE EQUATIONS
if j == 1
   GG(1:3,1) = -phi(1:3)';
end
%else
    k = k+(weight*E*area*jac)*(dphi(j,:)'*dphi(j,:));
    fbody = area*xg;
    f = f+(weight*fbody*jac)*phi(j,:)';
%end
end

% ENFORCE BOUNDARY CONDITIONS USING LAGRANGE MULTIPLIERS
q = [0];
m = [k GG;GG' zeros(1)];

% SOLVE FOR NODAL PARAMETERS
d = m\[f' q]';
u = d(1:nnodes);

result = phi*u;
for i=1:length(x)
    exact(i) =  0.5*(x(i)) - x(i)^3/6.0;
end
%%
plot(x,u);
title('valores u');
figure();
plot(x,result);
title('valores u*phi');
figure();
plot(x,exact);
title('valores exatos');
