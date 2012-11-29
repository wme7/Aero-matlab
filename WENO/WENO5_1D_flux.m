function [uun,uup] = WENO5_1D_flux(uu)
% 1D WENO3 reconstruction for 
% u_t + f_x = 0 where f(u) = a*u(x) 
% we assume "a" is or is not constant along x

%% NOTATION
% uu(1:9) nine points Flux values
%
% pjk  (j,k = 0,1,2,3,4) : the coefficients of the k-th degree term of the
% polynomial pj by weno5, e.g.,p12 represents the 2nd-degree term (x^2) of
% the polynomial p1(x) by weno5 
%
% dpjk (j,k = 0,1,2,3,4) : the coefficient of the k-th degree term of the
% derivative of the polynomial pj by weno5
%
% bj  (j = 0,1,2,3,4) : the smoothness indicator for the polynomial pj
%
% rnk (k = 0,1,2,3,4) : the optimal coefficients for the polynomials obtained via left stencils 
% rpk (k = 0,1,2,3,4) : the optimal coefficients for the polynomials obtained via right stencils
%
% wnk (k = 0,1,2,3,4) : the weights for the polynomials obtained via left stencils
% wpk (k = 0,1,2,3,4) : the weights for the polynomials obtained via right stencils

% Constants
dn0  = 1./126.; dn1  = 10./63.; dn2  = 10./21.; dn3  = 20./63.; dn4  = 5./126.;
dp0  = 5./126.; dp1  = 20./63.; dp2  = 10./21.; dp3  = 10./63.; dp4  = 1./126.;
epsilon  = 1E-6;

% Center value:
i = 5;

p00 = (- 71.*uu(i-4)+ 364.*uu(i-3)- 746.*uu(i-2)+ 684.*uu(i-1)+1689.*uu(i))/1920.;
p01 = (   9.*uu(i-4)-  50.*uu(i-3)+ 120.*uu(i-2)- 174.*uu(i-1)+  95.*uu(i))/48.;
p02 = (   7.*uu(i-4)-  36.*uu(i-3)+  74.*uu(i-2)-  68.*uu(i-1)+  23.*uu(i))/16.;
p03 = (   3.*uu(i-4)-  14.*uu(i-3)+  24.*uu(i-2)-  18.*uu(i-1)+   5.*uu(i))/12.;
p04 = (      uu(i-4)-   4.*uu(i-3)+   6.*uu(i-2)-   4.*uu(i-1)+      uu(i))/24.;
%
p10 = (   9.*uu(i-3)-  36.*uu(i-2)-  26.*uu(i-1)+2044.*uu(i)-  71.*uu(i+1))/1920.;
p11 = (-  5.*uu(i-3)+  30.*uu(i-2)-  84.*uu(i-1)+  50.*uu(i)+   9.*uu(i+1))/48.;
p12 = (-     uu(i-3)+   4.*uu(i-2)+   2.*uu(i-1)-  12.*uu(i)+   7.*uu(i+1))/16.;
p13 = (      uu(i-3)-   6.*uu(i-2)+  12.*uu(i-1)-  10.*uu(i)+   3.*uu(i+1))/12.;
p14 = (      uu(i-3)-   4.*uu(i-2)+   6.*uu(i-1)-   4.*uu(i)+      uu(i+1))/24.;
%
p20 = (   9.*uu(i-2)- 116.*uu(i-1)+2134.*uu(i)- 116.*uu(i+1)+   9.*uu(i+2))/1920.;
p21 = (   5.*uu(i-2)-  34.*uu(i-1)            +  34.*uu(i+1)-   5.*uu(i+2))/48.;
p22 = (-     uu(i-2)+  12.*uu(i-1)-  22.*uu(i)+  12.*uu(i+1)-      uu(i+2))/16.;
p23 = (-     uu(i-2)+   2.*uu(i-1)            -   2.*uu(i+1)+      uu(i+2))/12.;
p24 = (      uu(i-2)-   4.*uu(i-1)+   6.*uu(i)-   4.*uu(i+1)+      uu(i+2))/24.;
%
p30 = (- 71.*uu(i-1)+2044.*uu(i)-  26.*uu(i+1)-  36.*uu(i+2)+   9.*uu(i+3))/1920.;
p31 = (-  9.*uu(i-1)-  50.*uu(i)+  84.*uu(i+1)-  30.*uu(i+2)+   5.*uu(i+3))/48.;
p32 = (   7.*uu(i-1)-  12.*uu(i)+   2.*uu(i+1)+   4.*uu(i+2)-      uu(i+3))/16.;
p33 = (-  3.*uu(i-1)+  10.*uu(i)-  12.*uu(i+1)+   6.*uu(i+2)-      uu(i+3))/12.;
p34 = (      uu(i-1)-   4.*uu(i)+   6.*uu(i+1)-   4.*uu(i+2)+      uu(i+3))/24.;
%
p40 = (1689.*uu(i)+ 684.*uu(i+1)- 746.*uu(i+2)+ 364.*uu(i+3)-  71.*uu(i+4))/1920.;
p41 = (- 95.*uu(i)+ 174.*uu(i+1)- 120.*uu(i+2)+  50.*uu(i+3)-   9.*uu(i+4))/48.;
p42 = (  23.*uu(i)-  68.*uu(i+1)+  74.*uu(i+2)-  36.*uu(i+3)+   7.*uu(i+4))/16.;
p43 = (-  5.*uu(i)+  18.*uu(i+1)-  24.*uu(i+2)+  14.*uu(i+3)-   3.*uu(i+4))/12.;
p44 = (      uu(i)-   4.*uu(i+1)+   6.*uu(i+2)-   4.*uu(i+3)+      uu(i+4))/24.;
%
dp01  = ( 4.*p04)^2/448. + (3.*p03)^2/80. + (2.*p02)^2/12. + p01^2 +  p04*p02/5. + p03*p01/2.;
dp02  = (12.*p04)^2/80.  + (6.*p03)^2/12. + (2.*p02)^2              +  p04*p02*4.;
dp03  = (24.*p04)^2/12.  + (6.*p03)^2;
dp04  = (24.*p04)^2;
%
dp11  = ( 4.*p14)^2/448. + (3.*p13)^2/80. + (2.*p12)^2/12. + p11^2 +  p14*p12/5. + p13*p11/2.;
dp12  = (12.*p14)^2/80.  + (6.*p13)^2/12. + (2.*p12)^2              +  p14*p12*4.;
dp13  = (24.*p14)^2/12.  + (6.*p13)^2;
dp14  = (24.*p14)^2;
%
dp21  = ( 4.*p24)^2/448. + (3.*p23)^2/80. + (2.*p22)^2/12. + p21^2 +  p24*p22/5. + p23*p21/2.;
dp22  = (12.*p24)^2/80.  + (6.*p23)^2/12. + (2.*p22)^2              +  p24*p22*4.;
dp23  = (24.*p24)^2/12.  + (6.*p23)^2;
dp24  = (24.*p24)^2;
%
dp31  = ( 4.*p34)^2/448. + (3.*p33)^2/80. + (2.*p32)^2/12. + p31^2 +  p34*p32/5. + p33*p31/2.;
dp32  = (12.*p34)^2/80.  + (6.*p33)^2/12. + (2.*p32)^2              +  p34*p32*4.;
dp33  = (24.*p34)^2/12.  + (6.*p33)^2;
dp34  = (24.*p34)^2;
%
dp41  = ( 4.*p44)^2/448. + (3.*p43)^2/80. + (2.*p42)^2/12. + p41^2 +  p44*p42/5. + p43*p41/2.;
dp42  = (12.*p44)^2/80.  + (6.*p43)^2/12. + (2.*p42)^2              +  p44*p42*4.;
dp43  = (24.*p44)^2/12.  + (6.*p43)^2;
dp44  = (24.*p44)^2;
%
b0   = dp01 + dp02 + dp03 + dp04;
b1   = dp11 + dp12 + dp13 + dp14;
b2   = dp21 + dp22 + dp23 + dp24;
b3   = dp31 + dp32 + dp33 + dp34;
b4   = dp41 + dp42 + dp43 + dp44;
%
wn0  = dn0 / (b0 + epsilon)^3;
wn1  = dn1 / (b1 + epsilon)^3;
wn2  = dn2 / (b2 + epsilon)^3;
wn3  = dn3 / (b3 + epsilon)^3;
wn4  = dn4 / (b4 + epsilon)^3;
wsum = wn0 + wn1 + wn2 + wn3 + wn4; %+ epsilon;
wn0  = wn0 / wsum;
wn1  = wn1 / wsum;
wn2  = wn2 / wsum;
wn3  = wn3 / wsum;
wn4  = wn4 / wsum;

wp0  = dp0 / (b0 + epsilon)^3;
wp1  = dp1 / (b1 + epsilon)^3;
wp2  = dp2 / (b2 + epsilon)^3;
wp3  = dp3 / (b3 + epsilon)^3;
wp4  = dp4 / (b4 + epsilon)^3;
wsum = wp0 + wp1 + wp2 + wp3 + wp4; %+ epsilon;
wp0  = wp0 / wsum;
wp1  = wp1 / wsum;
wp2  = wp2 / wsum;
wp3  = wp3 / wsum;
wp4  = wp4 / wsum;
% Positive Flux u_i+1/2 (+)
dx   =-0.5;
uup  = wp0 * (p00 + p01*dx + p02*dx^2 + p03*dx^3 + p04*dx^4) ...
    + wp1 * (p10 + p11*dx + p12*dx^2 + p13*dx^3 + p14*dx^4) ...
    + wp2 * (p20 + p21*dx + p22*dx^2 + p23*dx^3 + p24*dx^4) ...
    + wp3 * (p30 + p31*dx + p32*dx^2 + p33*dx^3 + p34*dx^4) ...
    + wp4 * (p40 + p41*dx + p42*dx^2 + p43*dx^3 + p44*dx^4);
% Negative Flux u_i+1/2 (-)
dx   = 0.5;
uun  = wn0 * (p00 + p01*dx + p02*dx^2 + p03*dx^3 + p04*dx^4) ...
    + wn1 * (p10 + p11*dx + p12*dx^2 + p13*dx^3 + p14*dx^4) ...
    + wn2 * (p20 + p21*dx + p22*dx^2 + p23*dx^3 + p24*dx^4) ...
    + wn3 * (p30 + p31*dx + p32*dx^2 + p33*dx^3 + p34*dx^4) ...
    + wn4 * (p40 + p41*dx + p42*dx^2 + p43*dx^3 + p44*dx^4);
end