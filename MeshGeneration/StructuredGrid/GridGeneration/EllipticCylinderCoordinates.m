function [x y] = EllipticCylinderCoordinates(Xi,Eta,a)

% a is radius of ellips
% a > 0
if nargin == 2
    a = 1. ;
end

r = 1+Xi ;
s = pi*Eta ;

x = a*cosh(r)*cos(s) ;
y = a*sinh(r)*sin(s) ;