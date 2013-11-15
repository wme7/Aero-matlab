function a=Anonlin(u,i)
% A nonlinear non-decreasing function 
  if nargin<2 || (i==0),
     a=(u+0.25).*(u<=-0.25)+(u-0.25).*(u>=0.25);  %a=((u-0.5).^3)/(3.0);
  else
     a=1*(abs(u)>=0.25);%a=(u-0.5).^2;
  end