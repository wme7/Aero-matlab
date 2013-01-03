function res = hybrid(p,type)
% 	Used to switch from upwind to central scheme
n = length(p); res1 = ones(n,1);
if type == 1
  res1 = ones(n,1);
elseif type == 2
  res1 = zeros(n,1);
else
  for j = 1:n
    if abs(p(j)) < 1.9
      res1(j) = 0;
    elseif abs(p(j)) > 2
      res1(j) =  1;
    else
      res1(j) = 10*(abs(p(j))-1.9);
    end
  end
end
res = res1;
