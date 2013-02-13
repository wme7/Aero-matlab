function h = WENOtest2(u)
% Polynomial coef, From table:

% $c_{rj}$ for $x_{i+1/2}^{-}$
C00 =  1/3; C01 =  5/6; C02 = -1/6;
C10 = -1/6; C11 =  5/6; C12 =  1/3;
C20 =  1/3; C21 = -7/6; C22 = 11/6;

% % $c^{~}_{rj}$ = c_{r-1,j}$ for $x_{i+1/2}^{+}$
% c00 = 11/6; c01 = -7/6; c02 =  1/3;
% c10 =  1/3; c11 =  5/6; c12 = -1/6;
% c20 = -1/6; c21 =  5/6; c22 =  1/3;

% Polynomials p_r^{-} & p_r^{+}
%for i = 4:nx-3
i = 3; % reference of stencil

% coord: i - r +j
p0n = C00*u( i ) + C01*u(i+1) + C02*u(i+2);
p1n = C10*u(i-1) + C11*u( i ) + C12*u(i+1);
p2n = C20*u(i-2) + C21*u(i-1) + C22*u( i );

% Beta values
B0n = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(3*u( i )-4*u(i+1)+u(i+2))^2;
B1n = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-u(i+1))^2;
B2n = 13/12*(u(i-2)-2*u(i-1)+u( i ))^2 + 1/4*(u(i-2)-4*u(i-1)+3*u( i ))^2;

% Constants
d0n = 3/10; d1n = 6/10; d2n = 1/10; 
%d0p = 1/10; d1p = 6/10; d2p = 3/10;
epsilon = 1E-7;

% Non-normalized stencil weights 
alpha0n = d0n /(epsilon + B0n)^2;
alpha1n = d1n /(epsilon + B1n)^2;
alpha2n = d2n /(epsilon + B2n)^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% n-Weigths
w0n = alpha0n / alphasumn;
w1n = alpha1n / alphasumn;
w2n = alpha2n / alphasumn;



% % coord: (i+1)-((r-1)+1) + j
% p0p = c00*u(i+1) + c01*u(i+2) + c02*u(i+3);
% p1p = c10*u( i ) + c11*u(i+1) + c12*u(i+2);
% p2p = c20*u(i-1) + c21*u( i ) + c22*u(i+1);
% 
% % Beta values
% B0p = 13/12*(u(i+1)-2*u(i+2)+u(i+3))^2 + 1/4*(3*u(i+1)-4*u(i+2)+u(i+3))^2;
% B1p = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(u( i )-u(i+2))^2;
% B2p = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-4*u( i )+3*u(i+1))^2;
% 
% % Non-normalized stencil weights 
% alpha0p = d0p /(epsilon + B0p)^2;
% alpha1p = d1p /(epsilon + B1p)^2;
% alpha2p = d2p /(epsilon + B2p)^2;
% alphasump = alpha0p + alpha1p + alpha2p;
% 
% % p-Weigths
% w0p = alpha0p / alphasump;
% w1p = alpha1p / alphasump;
% w2p = alpha2p / alphasump;

%end

% un = u_{x+1/2}^{-} & up = u_{x+1/2}^{+}
un = w0n*p0n + w1n*p1n + w2n*p2n;
% up = w0p*p0p + w1p*p1p + w2p*p2p;

% Lax-Friedrichs, h(u): our final approximation of u_(i+1/2)
h = un; %0.5*(un+up-a*(un-up));
