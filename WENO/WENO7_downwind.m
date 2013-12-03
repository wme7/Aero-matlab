function h = WENO7_downwind(u)
% *************************************************************************
% This subroutine 'downwinds' our information from right to left.
% This is done by assuming u(i) are the cell's averages flux information in
% our domain (see domain ref), the cell boundary fluxes are computed using
% WENO3 (k=2) reconstruction. 
%
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: h(i) = $u_{i-1/2}^{+}$
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
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]

% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 

% Polynomials p_r^{+} 
 i = 4; % is the central reference for the following stencil, S{i};

% Polynomial coefcients (Taken from table)
% $\tilde c_{rj}$ = c_{r-1,j}$ for $x_{i+1/2}^{+}$
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
up = w0p*p0p + w1p*p1p + w2p*p2p + w3p*p3p;

% Numerical Flux
% We assume negative advection speed,(a < 0),therefore wind blows from
% right to left. We would use u_{i-1/2}^{+}, (un), for the numerical flux;
h = up;
