function [x y] = BipolarCoordinates(Xi,Eta,a) 

% a is radii that defies height of the domain
% a > 0

% Check the input arguments and set default values

if nargin == 2
    a = 1. ;
end
r = Xi ;
s = pi*(Eta-1/2) ;

x = (a*sinh(r))/(cosh(r)+cos(s)) ;
y = (a*sin(s))/(cosh(r)+cos(s)) ;