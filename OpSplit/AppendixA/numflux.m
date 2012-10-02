function F=numflux(flux,u,method,lambda)
%
% Returns the numerical flux function used in a finite volume
% method for the conservation law with flux "flux". 
%
% method is either Lax-Friedrichs ('lxf'), Lax-Wendroff ('lxw'), MacCormac
% ('mcc') or upwind ('upwind'). Other methods can easily be added
%
S=size(u);
n=S(1)-1;
f=feval(flux,u);
switch lower(method)
	case('lxf')
		F=(u(1:n,:)-u(2:n+1,:))/(2*lambda)+0.5*(f(1:n,:)+f(2:n+1,:));
	case('lxw')
		uh=0.5*(u(1:n,:)+u(2:n+1,:)-lambda*(f(2:n+1,:)-f(1:n,:)));
		F=feval(flux,uh);
	case('mcc')
		uh=u(2:n+1,:)-lambda*(f(2:n+1,:)-f(1:n,:));
		F=0.5*(feval(flux,uh)+f(1:n,:));
	case('upwind')
		F=f(1:n,:);
	otherwise
		mess=strcat(method,' not implemented in numflux.m');
		error(mess);
end;