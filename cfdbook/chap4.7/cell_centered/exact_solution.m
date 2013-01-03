function ye = exact_solution(x,y,D)
z = (1/sqrt(2-x))*(exp(-y*y/(4*D*(2-x))) + exp(-(2-y)*(2-y)/(4*D*(2-x))));
ye = z;
