function [x y] = ParabolicCylinderCoordinates(Xi,Eta)

r = 1+Xi ;
s = 1+Eta ;

x = 1/2*(r^2-s^2) ;
y = r*s ;