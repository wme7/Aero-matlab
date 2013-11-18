% use umP'th order interpolation
umP = 3;

% load triangle nodes for umP+1'th order interpolation
load(sprintf('um_n%d.mat', umP+1))

% use symbolic computation 
r = sym('r');
s = sym('s');

% symbolically compute the polynomial expansions
% for the Proriol-Koornwinder-Dubiner-Owens
% basis up to umP'th total order
sk = 1;
for m=0:umP
  for n=0:m
    pkdo(sk) = umSYMBPKDO2d(n,m-n);
    sk=sk+1;
  end
end
Nmembers = sk-1;

% symbolically compute PKDO derivative matrix
for n=1:Nmembers
  dpkdodr(n,1) = diff(pkdo(n), r);
  dpkdods(n,1) = diff(pkdo(n), s);
end
dig = 32;
digits(dig)

if(0)
% symbollically compute PKDO 1D mass matrices for r and s integration
for n=1:Nmembers
  for m=1:Nmembers
    integrand = pkdo(n)*pkdo(m);

    % integrate product of two PKDO members
    % over triangle and simplify
    facemm1(n,m) = simple(int(integrand, r, -1, 1));
    facemm2(n,m) = simple(int(integrand, s, -1, 1));
  end
end

% evaluate face mass matrices

s = vpa('-1', dig); fmm1  = eval(facemm1);
r = vpa( '1', dig);  fmm2  = eval(facemm2);
r = vpa('-1', dig);  fmm3  = eval(facemm2);

% compute Nodal Face Mass Matrices 
FaceMassMatrix1  = double(transpose(vdm)\fmm1/vdm);
FaceMassMatrix2  = double(transpose(vdm)\fmm2/vdm);
FaceMassMatrix3  = double(transpose(vdm)\fmm3/vdm);

end

vdm      = sym(zeros(umNpts, Nmembers), 'd');
dvdmdr = sym(zeros(umNpts, Nmembers), 'd');
dvdmds = sym(zeros(umNpts, Nmembers), 'd');

% Build generalized Vandermonde matrices (and derivatives)
for n=1:umNpts
  r = vpa(sym(umR(n)), dig);
  s = vpa(sym(umS(n)), dig);
  for m=1:Nmembers
    vdm(n,m)      = sym(eval(pkdo(m)   ));
    dvdmdr(n,m) = sym(eval(dpkdodr(m)));
    dvdmds(n,m) = sym(eval(dpkdods(m)));
  end
end

% compute Nodal Mass Matrix
MassMatrix = double(inv(transpose(vdm))/vdm);


% compute Nodal derivative matrix
Dr = double(dvdmdr/vdm);
Ds = double(dvdmds/vdm);
