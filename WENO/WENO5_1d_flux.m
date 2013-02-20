function [hn,hp] = WENO5_1d_flux(v,u)
% Compute numerical fluxes at cell 'i' interfaces.
% Input:  v(1:7) = positive fluxes - cells average values
%         u(1:7) = negative fluxes - cells average values
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
 i = 4; % is the central for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $c_{rj}$ for $u_{i+1/2}^{-}$
% C00 =  1/4; C01 =13/12; C02 =-5/12; C03 = 1/12;
% C10 =-1/12; C11 = 7/12; C12 = 7/12; C13 =-1/12;
% C20 = 1/12; C21 = -5/6; C22 =13/12; C23 = 1/4 ;
% C30 = -1/4; C31 =13/12; C32=-23/12; C33 =25/12;

C00 = -1/4; C31 =13/12; C32=-23/12; C33 =25/12;
C10 = 1/12; C21 = -5/6; C22 =13/12; C23 = 1/4 ;
C20 =-1/12; C11 = 7/12; C12 = 7/12; C13 =-1/12;
C30 =  1/4; C01 =13/12; C02 =-5/12; C03 = 1/12;

% coord: i - r +j
p0n = C00*v(i-3) + C01*v(i-2) + C02*v(i-1) + C03*v( i );
p1n = C10*v(i-2) + C11*v(i-1) + C12*v( i ) + C13*v(i+1);
p2n = C20*v(i-1) + C21*v( i ) + C22*v(i+1) + C23*v(i+2);
p3n = C30*v( i ) + C31*v(i+1) + C32*v(i+2) + C33*v(i+3);

% Smooth Indicators, Beta factors
B0n = v(i-3)*(  547*v(i-3)- 3882*v(i-2)+ 4642*v(i-1) - 1854*v( i )) ...
    + v(i-2)*( 7043*v(i-2)-17246*v(i-1)+ 7042*v( i )) ...
    + v(i-1)*(11003*v(i-1)- 9402*v( i )) ...
    + 2107*v( i )^2;
B1n = v(i-2)*(  267*v(i-2)- 1642*v(i-1)+ 1602*v( i ) -  494*v(i+1)) ...
    + v(i-1)*( 2843*v(i-1)- 5966*v( i )+ 1922*v(i+1)) ...
    + v( i )*( 3443*v( i )- 2522*v(i+1)) ...
    + 547*v(i+1)^2;
B2n = v(i-1)*(  547*v(i-1)- 2522*v( i )+ 1922*v(i+1) -  494*v(i+2)) ...
    + v( i )*( 3443*v( i )- 5966*v(i+1)+ 1602*v(i+2)) ...
    + v(i+1)*( 2843*v(i+1)- 1642*v(i+2)) ...
    + 267*v(i+2)^2;
B3n = v( i )*( 2107*v( i )- 9402*v(i+1)+ 7042*v(i+2) - 1854*v(i+3)) ...
    + v(i+1)*(11003*v(i+1)-17246*v(i+2)+ 4642*v(i+3)) ...
    + v(i+2)*( 7043*v(i+2)- 3882*v(i+3)) ...
    + 547*v(i+3)^2;

% Constants
d0n = 1/35; d1n = 12/35; d2n = 18/35; d3n = 4/35;
epsilon = 1E-7;

% Alpha weights 
alpha0n = d0n /(epsilon + B0n)^2;
alpha1n = d1n /(epsilon + B1n)^2;
alpha2n = d2n /(epsilon + B2n)^2;
alpha3n = d3n /(epsilon + B3n)^2;
alphasumn = alpha0n + alpha1n + alpha2n + alpha3n;

% ENO stencils weigths
w0n = alpha0n / alphasumn;
w1n = alpha1n / alphasumn;
w2n = alpha2n / alphasumn;
w3n = alpha2n / alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n*p0n + w1n*p1n + w2n*p2n + w3n*p3n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 

% Polynomials p_r^{+} 
 i = 4; % is the central reference for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $\tilde c_{rj}$ = c_{r-1,j}$ for $x_{i+1/2}^{+}$
% c00 =25/12; c01=-23/12; c02 =13/12; c03 = -1/4; 
% c10 =  1/4; c11 =13/12; c12 =-5/12; c13 = 1/12;
% c20 =-1/12; c21 = 7/12; c22 = 7/12; c23 =-1/12;
% c30 = 1/12; c31 = -5/6; c32 =13/12; c33 =  1/4;

c00 = 1/12; c31 = -5/6; c32 =13/12; c33 =  1/4;
c10 =-1/12; c21 = 7/12; c22 = 7/12; c23 =-1/12;
c20 =  1/4; c11 =13/12; c12 =-5/12; c13 = 1/12;
c30 =25/12; c01=-23/12; c02 =13/12; c03 = -1/4;

% coord: (i)-((r-1)+1) + j
p0p = c00*u(i-3) + c01*u(i-2) + c02*u(i-1) + c03*u( i );
p1p = c10*u(i-2) + c11*u(i-1) + c12*u( i ) + c13*u(i+1);
p2p = c20*u(i-1) + c21*u( i ) + c22*u(i+1) + c23*u(i+2);
p3p = c30*u( i ) + c31*u(i+1) + c32*u(i+2) + c33*u(i+3);

% Smooth Indicators, Beta factors
B0p = u(i-3)*(  547*u(i-3)- 3882*u(i-2)+ 4642*u(i-1) - 1854*u( i )) ...
    + u(i-2)*( 7043*u(i-2)-17246*u(i-1)+ 7042*u( i )) ...
    + u(i-1)*(11003*u(i-1)- 9402*u( i )) ...
    + 2107*u( i )^2;
B1p = u(i-2)*(  267*u(i-2)- 1642*u(i-1)+ 1602*u( i ) - 494*u(i+1)) ...
    + u(i-1)*( 2843*u(i-1)- 5966*u( i )+ 1922*u(i+1)) ...
    + u( i )*( 3443*u( i )- 2522*u(i+1)) ...
    + 547*u(i+1)^2;
B2p = u(i-1)*(  547*u(i-1)- 2522*u( i )+ 1922*u(i+1) - 494*u(i+2)) ...
    + u( i )*( 3443*u( i )- 5966*u(i+1)+ 1602*u(i+2)) ...
    + u(i+1)*( 2843*u(i+1)- 1642*u(i+2)) ...
    + 267*u(i+2)^2;
B3p = u( i )*( 2107*u( i )- 9402*u(i+1)+ 7042*u(i+2) - 1854*u(i+3)) ...
    + u(i+1)*(11003*u(i+1)-17246*u(i+2)+ 4642*u(i+3)) ...
    + u(i+2)*( 7043*u(i+2)- 3882*u(i+3)) ...
    + 547*u(i+3)^2;

% Constants
d0p = 1/35; d1p = 12/35; d2p = 18/35; d3p = 4/35;
epsilon = 1E-7;

% Alpha weights 
alpha0p = d0p /(epsilon + B0p)^2;
alpha1p = d1p /(epsilon + B1p)^2;
alpha2p = d2p /(epsilon + B2p)^2;
alpha3p = d3p /(epsilon + B3p)^2;
alphasump = alpha0p + alpha1p + alpha2p + alpha3p;

% ENO stencils weigths
w0p = alpha0p / alphasump;
w1p = alpha1p / alphasump;
w2p = alpha2p / alphasump;
w3p = alpha3p / alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p*p0p + w1p*p1p + w2p*p2p + w3p*p3p;

% Total flux
%h = hn + hp;