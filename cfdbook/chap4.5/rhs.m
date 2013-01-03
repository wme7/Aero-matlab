function res = rhs(y,a,b,alpha,pe)
% Right-hand side of differential equation
res = zeros(length(y),1); ppi = pi; bet = 1/(exp(pe)-1); n = length(y);
for j=1:n
  uu = 1 - b*sin(ppi*y(j));
  vi = a*(sin(alpha*ppi*y(j)) - sin(alpha*ppi)) + bet*(exp(y(j)*pe) - 1);
  ress = uu*(a*alpha*ppi*cos(alpha*ppi*y(j)) + bet*pe*exp(pe*y(j)));
  ress = ress - vi*b*ppi*cos(ppi*y(j));
  ress = ress + (a/pe)*(alpha*ppi)*(alpha*ppi)*sin(alpha*ppi*y(j));
  ress = ress - bet*pe*exp(pe*y(j));
  res(j) = ress;
end
