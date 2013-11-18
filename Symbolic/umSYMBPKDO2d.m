
function pkdo = umSYMBPKDO2d(i,j)

  r = sym('r');
  s = sym('s');
  
  a = 2*(1+r)/(1-s) - 1;
  b = s; 
  
  pkdo = umSYMBJACOBI1d(a,i,0,0);
  pkdo = pkdo.*umSYMBJACOBI1d(b,j,2*i+1,0).*((0.5*(1-b)).^i);

  % normalize (sets L2 norm of each basis function to be 1)
  pkdo = pkdo*sqrt((i+0.5)*(i+j+1));
  
  % simplify
  pkdo = simple(pkdo);  

