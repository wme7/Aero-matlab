function [hn,hp] = WENO3_1d_flux(v,u)
% Compute numerical fluxes at cell 'i' interfaces.
% Input:  v(1:6) = positive fluxes - cells average values
%         u(1:6) = negative fluxes - cells average values
% Output: un(1)  = Numerical flux @ x_{i+1/2}^(-) | right flux
%         up(1)  = Numerical flux @ x_{i-1/2}^(+) | left flux
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |           |h-       h+|           |
%                          i-1/2       i+1/2

%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$

% Polynomials p_r^{-} 
 i = 3; % is the central for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $c_{rj}$ for $u_{i+1/2}^{-}$
C00 =  1/3; C01 = -7/6; C02 = 11/6;
C10 = -1/6; C11 =  5/6; C12 =  1/3;
C20 =  1/3; C21 =  5/6; C22 = -1/6;

% coord: i - r +j
p0n = C00*v(i-2) + C01*v(i-1) + C02*v( i );
p1n = C10*v(i-1) + C11*v( i ) + C12*v(i+1);
p2n = C20*v( i ) + C21*v(i+1) + C22*v(i+2);

% Smooth Indicators, Beta factors
B0n = 13/12*(v(i-2)-2*v(i-1)+v( i ))^2 + 1/4*(v(i-2)-4*v(i-1)+3*v( i ))^2;
B1n = 13/12*(v(i-1)-2*v( i )+v(i+1))^2 + 1/4*(v(i-1)-v(i+1))^2;
B2n = 13/12*(v( i )-2*v(i+1)+v(i+2))^2 + 1/4*(3*v( i )-4*v(i+1)+v(i+2))^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; 
epsilon = 1E-7;

% Alpha weights 
alpha0n = d0n /(epsilon + B0n)^2;
alpha1n = d1n /(epsilon + B1n)^2;
alpha2n = d2n /(epsilon + B2n)^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n / alphasumn;
w1n = alpha1n / alphasumn;
w2n = alpha2n / alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n*p0n + w1n*p1n + w2n*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 

% Polynomials p_r^{+} 
 i = 3; % is the central reference for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $\tilde c_{rj}$ = c_{r-1,j}$ for $x_{i+1/2}^{+}$
c00 = -1/6; c01 =  5/6; c02 =  1/3;
c10 =  1/3; c11 =  5/6; c12 = -1/6;
c20 = 11/6; c21 = -7/6; c22 =  1/3;

% coord: (i)-((r-1)+1) + j
p0p = c00*u(i-2) + c01*u(i-1) + c02*u( i );
p1p = c10*u(i-1) + c11*u( i ) + c12*u(i+1);
p2p = c20*u( i ) + c21*u(i+1) + c22*u(i+2);

% Smooth Indicators, Beta factors
B0p = 13/12*(u(i-2)-2*u(i-1)+u( i ))^2 + 1/4*(u(i-2)-4*u(i-1)+3*u(i))^2;
B1p = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-u(i+1))^2;
B2p = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(3*u(i)-4*u(i+1)+u(i+2))^2;

% Constants
d0p = 1/10; d1p = 6/10; d2p = 3/10;
epsilon = 1E-7;

% Alpha weights 
alpha0p = d0p /(epsilon + B0p)^2;
alpha1p = d1p /(epsilon + B1p)^2;
alpha2p = d2p /(epsilon + B2p)^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p / alphasump;
w1p = alpha1p / alphasump;
w2p = alpha2p / alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p*p0p + w1p*p1p + w2p*p2p;

% Total flux
%h = hn + hp;