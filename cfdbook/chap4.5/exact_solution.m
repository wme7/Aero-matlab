function res = exact_solution(y,a,alpha,pe)
%       Exact solution: v = a{sin(alpha.pi.x) - sin(alpha.pi)} +
%       \frac{exp[(x-1)Pe]-exp(-Pe)]}{1-exp(-Pe)}
n = length(y); res1 = zeros(n,1);
for j = 1:n
  res1(j) = a*(sin(alpha*pi*y(j)) - sin(alpha*pi)) + ...
  (exp((y(j)-1)*pe)-exp(-pe))/(1-exp(-pe));
end
res = res1;
