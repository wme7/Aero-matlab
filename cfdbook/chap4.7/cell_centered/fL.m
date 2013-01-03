function ye = fL(y, D, hom)
% Neumann condition at outflow
z =  2^(-2.5)*(y*y*0.25/D - 1)* exp(-y*y/(8*D));
z = z + 2^(-2.5)*((2-y)*(2-y)*0.25/D - 1)*exp(-(2-y)*(2-y)/(8*D));
ye = hom*z;
