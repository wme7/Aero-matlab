function h = WENO7_upwind(v)
% *************************************************************************
% This subroutine 'upwinds' our information from left to right.
% This is done by assuming u(i) are the cell's averages flux information in
% our domain (see domain ref), the cell boundary fluxes are computed using
% WENO3 (k=2) reconstruction. 
%
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: h(i) = $u_{i+1/2}^{-}$
%
% Based on:
% C.W. Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
% *************************************************************************
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |___________S3__________|
%                               |                       |
%                       |___________S2__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
% WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]

% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$

% Polynomials p_r^{-} 
 i = 4; % is the central for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $c_{rj}$ for $u_{i+1/2}^{-}$
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
vn = w0n*p0n + w1n*p1n + w2n*p2n + w3n*p3n;

% Numerical Flux
% We assume advection speed positive (a > 0) therefore wind blows from left
% to right. We would use u_{i+1/2}^{-}, (un), for the numerical flux;
h = vn;
