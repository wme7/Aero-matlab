% .............................................................     
    function [shape,naturalDerivatives]=shapeFunctionQ12(xi,eta)

    % shape function and derivatives for Q12 elements
    % shape : Shape functions
    % naturalDerivatives: derivatives w.r.t. xi and eta 
    % xi, eta: natural coordinates (-1 ... +1)
   
    shape=[1/32*(1-xi)*(1-eta)*(-10+9*(xi^2+eta^2));
           1/32*(1+xi)*(1-eta)*(-10+9*(xi^2+eta^2));
           1/32*(1+xi)*(1+eta)*(-10+9*(xi^2+eta^2));
           1/32*(1-xi)*(1+eta)*(-10+9*(xi^2+eta^2));
           9/32*(1-xi^2)*(1-1/3*xi)*(1-eta);
           9/32*(1-xi^2)*(1+1/3*xi)*(1-eta);
           9/32*(1+xi)*(1-eta^2)*(1-1/3*eta);
           9/32*(1+xi)*(1-eta^2)*(1+1/3*eta);
           9/32*(1-xi^2)*(1+1/3*xi)*(1+eta);
           9/32*(1-xi^2)*(1-1/3*xi)*(1+eta);
           9/32*(1-xi)*(1-eta^2)*(1+1/3*eta);
           9/32*(1-xi)*(1-eta^2)*(1-1/3*eta)];
                   
                       
       naturalDerivatives=...
       [9/16*(1-eta)*(1-xi)*xi-(1/32)*(1-eta)*(-10+9*(eta^2+xi^2)),9/16*(1-eta)*eta*(1-xi)-1/32*(1-xi)*(-10+9*(eta^2+xi^2));
       9/16*(1-eta)*xi*(1+xi)+(1/32)*(1-eta)*(-10+9*(eta^2+xi^2)),9/16*(1-eta)*eta*(1+xi)-(1/32)*(1+xi)*(-10+9*(eta^2+xi^2)) ;
       9/16*(1+eta)*xi*(1+xi)+(1/32)*(1+eta)*(-10+9*(eta^2+xi^2)),9/16*eta*(1+eta)*(1+xi)+(1/32)*(1+xi)*(-10+9*(eta^2+xi^2)) ; 
       9/16*(1+eta)*(1-xi)*xi-(1/32)*(1+eta)*(-10+9*(eta^2+xi^2)),9/16*eta*(1+eta)*(1-xi)+(1/32)*(1-xi)*(-10+9*(eta^2+xi^2)) ;
       -(9/16)*(1-eta)*(1-3*xi)*xi-(27/32)*(1-eta)*(1-xi^2),-(9/32)*(1-3*xi)*(1-xi^2);
       -(9/16)*(1-eta)*xi*(1+3*xi)+(27/32)*(1-eta)*(1-xi^2),-(9/32)*(1+3*xi)*(1-xi^2);
       (9/32)*(1-3*eta)*(1-eta^2),-(9/16)*(1-3*eta)*eta*(1+xi)-(27/32)*(1-eta^2)*(1+xi);
       (9/32)*(1+3*eta)*(1-eta^2),-(9/16)*eta*(1+3*eta)*(1+xi)+(27/32)*(1-eta^2)*(1+xi);
     -(9/16)*(1+eta)*xi*(1+3*xi)+(27/32)*(1+eta)*(1-xi^2),(9/32)*(1+3*xi)*(1-xi^2);
       -(9/16)*(1+eta)*(1-3*xi)*xi-(27/32)*(1+eta)*(1-xi^2),(9/32)*(1-3*xi)*(1-xi^2);
       -(9/32)*(1+3*eta)*(1-eta^2),-(9/16)*eta*(1+3*eta)*(1-xi)+(27/32)*(1-eta^2)*(1-xi);
       -(9/32)*(1-3*eta)*(1-eta^2),-(9/16)*(1-3*eta)*eta*(1-xi)-(27/32)*(1-eta^2)*(1-xi)];
            
    end % end function shapeFunctionQ12
  
