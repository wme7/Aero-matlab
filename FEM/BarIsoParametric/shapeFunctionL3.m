% .............................................................     
    function [shape,naturalDerivatives]=shapeFunctionL3(xi)

    % shape function and derivatives for L2 elements
    % shape : Shape functions
    % naturalDerivatives: derivatives w.r.t. xi 
    % xi: natural coordinates (-1 ... +1)
   
    shape = [xi*(xi-1), 2*(1-xi*xi), xi*(1+xi)]/2;
    
    naturalDerivatives = [2*xi-1, -4*xi, 2*xi+1]/2;
             
    end % end function shapeFunctionL2
  
