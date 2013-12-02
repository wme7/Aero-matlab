function h = WENO_downwind(u)
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

% Polynomial coefcients (Taken from table)
% $\tilde c_{rj}$ = c_{r-1,j}$ for $x_{i+1/2}^{+}$
c00 = -1/6; c01 =  5/6; c02 =  1/3;
c10 =  1/3; c11 =  5/6; c12 = -1/6;
c20 = 11/6; c21 = -7/6; c22 =  1/3;

% Polynomials p_r^{-} 
i = 3; % central reference of WENO stencil S{i}

% coord: (i)-((r-1)+1) + j
p0p = c00*u(i-2) + c01*u(i-1) + c02*u( i );
p1p = c10*u(i-1) + c11*u( i ) + c12*u(i+1);
p2p = c20*u( i ) + c21*u(i+1) + c22*u(i+2);

% Smooth Indicators, Beta factors
B0p = 13/12*(u(i-2)-2*u(i-1)+u( i ))^2 + 1/4*(u(i-2)-4*u(i-1)+3*u(i))^2;
B1p = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-u(i+1))^2;
B2p = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(3*u(i)-4*u(i+1)+u(i+2))^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10;
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
up = w0p*p0p + w1p*p1p + w2p*p2p;

% Numerical Flux
% We assume negative advection speed,(a < 0),therefore wind blows from
% right to left. We would use u_{i-1/2}^{+}, (un), for the numerical flux;
h = up;
