function u=initialfunc(x,k)
% initial function
if nargin<2,
	k=5;
end;
y=x-sign(x);
ekx=exp(-k*y);
u=1./(ekx+1);
