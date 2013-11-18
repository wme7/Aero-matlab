
r = sym('r');
s = sym('s');

% print out the polynomial expansions for the Proriol-Koornwinder-Dubiner-Owens
% basis up to total order 4 

sk = 1;
for m=0:4

  for n=0:m
    
    disp(sprintf('phi_{%d,%d} = ', n, m-n))
    res = umSYMBPKDO2d(n,m-n);

    pkdo(sk) = res; sk=sk+1;

  end
end

Nmembers = sk-1;

% demonstrate orthogonality

for n=1:Nmembers
  for m=1:Nmembers
    
    integrand = pkdo(n)*pkdo(m);

    % integrate product of two PKDO members over triangle
    massmatrix(n,m) = int(int(integrand, r, -1, -s), s, -1, 1);
  end
end

% display PKDO mass matrix
massmatrix


% compute PKDO derivative matrix

for n=1:Nmembers

  dpkdodr(n,1) = diff(pkdo(n), r);
  dpkdods(n,1) = diff(pkdo(n), s);
  
end

% set high number of digits for symbolic evaluations

load um_n10

for n=1:umNpts
  
  r = umR(n);
  s = umS(n);
  
  for m=1:Nmembers
    vdm(n,m)    = eval(pkdo(m)   );
    dvdmdr(n,m) = eval(dpkdodr(m));
    dvdmds(n,m) = eval(dpkdods(m));
  end
end

% compute Nodal Mass Matrix
MassMatrix = inv(transpose(vdm))/vdm;

% compute derivative matrix
Dr = dvdmdr/vdm;
Ds = dvdmds/vdm;
