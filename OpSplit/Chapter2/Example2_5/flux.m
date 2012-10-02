function f=flux(u,i)
if nargin<2 || (i==0)
   f=0.5*u.^2;    % Burgers' flux
elseif (i==1)
   f=u;           % Linear flux
else
	error(' Unknown option in flux ');
end;