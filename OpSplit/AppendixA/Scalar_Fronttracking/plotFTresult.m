function plotFTresult(u,x);
%
%  Plots the piecewise constant function u with discontinuities at x
%
[nc n]=size(u);   % nc = number of components
dp=0.8/nc;dpp=0.05;
for j=1:nc
  yhigh=max(u(j,:)); ylow=min(u(j,:));
  if ylow==yhigh,
    ylow=yhigh-0.5;
    yhigh=ylow+1;
  end;
  dy=0.05*(yhigh-ylow);
  if x(1)==x(n-1),
    xlow=x(1)-0.5;
    xhigh=xlow+1.0;
  else
    xhigh=x(n-1);
    xlow=x(1);
  end;
  dx=0.025*(xhigh-xlow);
  axis([xlow-dx,xhigh+dx,ylow-dy,yhigh+dy]);
  set(gca,'Box','on');
  x1=xlow-dx;
  for i=1:n-1
    line([x1 x(i)],[u(j,i) u(j,i)]);
    x1=x(i);
    line([x1 x1],[u(j,i) u(j,i+1)]);
  end;
  line([x1 xhigh+dx],[u(j,n),u(j,n)]);
end;
