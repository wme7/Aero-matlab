function h = WENO_upwind(u)
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
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]

% Polynomial coefcients (Taken from table)
% $c_{rj}$ for $u_{i+1/2}^{-}$
C00 =  1/3; C01 = -7/6; C02 = 11/6;
C10 = -1/6; C11 =  5/6; C12 =  1/3;
C20 =  1/3; C21 =  5/6; C22 = -1/6;

% Polynomials p_r^{-} 
i = 3; % central reference of WENO stencil S{i}

% coord: i - r +j
p0n = C00*u(i-2) + C01*u(i-1) + C02*u( i );
p1n = C10*u(i-1) + C11*u( i ) + C12*u(i+1);
p2n = C20*u( i ) + C21*u(i+1) + C22*u(i+2);

% Smooth Indicators, Beta factors
B0n = 13/12*(u(i-2)-2*u(i-1)+u( i ))^2 + 1/4*(u(i-2)-4*u(i-1)+3*u( i ))^2;
B1n = 13/12*(u(i-1)-2*u( i )+u(i+1))^2 + 1/4*(u(i-1)-u(i+1))^2;
B2n = 13/12*(u( i )-2*u(i+1)+u(i+2))^2 + 1/4*(3*u( i )-4*u(i+1)+u(i+2))^2;

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

% Compute aproximation $u_{i+1/2}^{-}$
un = w0n*p0n + w1n*p1n + w2n*p2n;

% Numerical Flux
% We assume advection speed positive (a > 0) therefore wind blows from left
% to right. We would use u_{i+1/2}^{-}, (un), for the numerical flux;
h = un;
