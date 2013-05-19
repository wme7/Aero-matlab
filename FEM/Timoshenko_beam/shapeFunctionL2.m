% .............................................................
function [shape,naturalDerivatives]=shapeFunctionL2(xi)
% shape function and derivatives for L2 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi
% xi: natural coordinates (-1 ... +1)
shape=([1-xi,1+xi]/2)';
naturalDerivatives=[-1;1]/2;
end % end function shapeFunctionL2