% .............................................................     
    function [shape,naturalDerivatives]=shapeFunctionQ4(xi,eta)

    % shape function and derivatives for Q4 elements
    % shape : Shape functions
    % naturalDerivatives: derivatives w.r.t. xi and eta 
    % xi, eta: natural coordinates (-1 ... +1)
   
    shape=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
            (1+xi)*(1+eta);(1-xi)*(1+eta)];
    naturalDerivatives=...
           1/4*[-(1-eta), -(1-xi);1-eta,    -(1+xi);
		         1+eta,      1+xi;-(1+eta),   1-xi];

    end % end function shapeFunctionQ4
  
