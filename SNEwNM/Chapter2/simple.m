function [f,jac]=simple(x)
% SIMPLE simple two-dimensional problem with interesting
%        global convergence behavior.
%        [F, JAC]=SIMPLE(X) are the function and Jacobian.
%
f=zeros(2,1);
f(1)=x(1)*x(1)+x(2)*x(2) - 2;
f(2)=exp(x(1)-1) + x(2)*x(2) - 2;
%
% Return the Jacobian if it's needed.
%
if nargout == 2
	jac=[2*x(1), 2*x(2); exp(x(1)-1), 2*x(2)];
end
