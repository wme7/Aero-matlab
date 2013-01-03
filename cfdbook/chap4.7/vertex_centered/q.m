function ye = q(x,y,D)
% right-hand-side function corresponding to exact solution
z = - 0.75*((2-x)^(-4.5)) * (D*(2-x)^2 - y*y*(2-x) + (y^4)/(12*D))*...
  exp(-y*y/(4*D*(2-x)));
z = z - 0.75*((2-x)^(-4.5)) * (D*(2-x)^2 - (2-y)*(2-y)*(2-x) +...
 ((2-y)^4)/(12*D)) * exp(-(2-y)*(2-y)/(4*D*(2-x)));
ye = z;
