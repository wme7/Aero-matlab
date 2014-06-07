% Computing states across shockwaves using 
% Rankine-Hugoniot (RH) condition
% Manuel Diaz, NTU, 2014.01.21

% Notation
% r: density
% u: velocity in x
% v: velocity in y
% p: presure

gamma = 5/3;

% state 1, ahead of the moving shock,
r1 = 1.4;
u1 = 0.0;
v1 = 0.0;
p1 = 1.0;

% Set incident Shock Mach number (Ms),
Ms = 10.0; 

% precompute Tau and c1,
c1 = sqrt(gamma*p1/r1); Tau = (gamma+1)/(gamma-1);

% state 2, behing of the shock,
p2 = p1*(2*gamma*Ms^2-(gamma-1))/(gamma+1);
r2 = r1*(Tau*(p2/p1) + 1)/(Tau + (p2/p1));
u2 = Ms*(1-((gamma-1)*Ms^2 + 2)/((gamma+1)*Ms^2))*c1;
v2 = 0.0;

disp(r2)
disp(u2)
disp(v2)
disp(p2)
