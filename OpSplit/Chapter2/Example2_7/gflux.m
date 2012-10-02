function f=gflux(u,i)
if nargin<2 || (i==0)
	% give the flux 
	f=-(1-u).*u;
elseif (i==1)
    % the derivative
	f=-(1-2*u);
else
	error(' Unknown option in gflux ');
end
% More advanced model of two-phase gravity flow
% un = u.^2;
% uw = (1-u).^2;
% if nargin<2 || (i==0)
% 	% give the flux 
% 	f=-un.*uw./(un + uw);
% elseif (i==1)
%     % the derivative 
% 	f=(-2*u.*uw.^2 + 2*un.^2.*(1-u))./(un+uw).^2;
% else
% 	error(' Unknown option in fflux ');
% end