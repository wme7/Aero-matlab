 % Copyright 2001, Brown University, Providence, Rhode Island.
 %
 % All Rights Reserved
 % 
 % Permission to use this software for noncommercial research and
 % educational purposes is hereby granted without fee.
 % Redistribution, sale, or incorporation of this software into a
 % commercial product is prohibited.
 % 
 % BROWN UNIVERSITY DISCLAIMS ANY AND ALL WARRANTIES WITH REGARD TO
 % THIS SOFTWARE,INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 % AND FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN
 % UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 % DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE,
 % DATA OR PROFITS.

% -----------------------------------------------------------------
%   jacobi() - jacobi polynomials 
%   
%   Get a vector 'poly' of values of the n_th order Jacobi polynomial
%   P^(alpha,beta)_n(z) alpha > -1, beta > -1 at the np points in z
%   -----------------------------------------------------------------

function jf = JACOBI1D(z, n, alpha, beta)

  dims = size(z);
  jf   = zeros(dims);
  one = 1.0;
  two = 2.0;

  if(n == 0)
   jf(:) = one;
  elseif (n == 1)
   jf(:) = 0.5*(alpha - beta + (alpha + beta + two)*z);
  else
   two = 2.0;
   apb = alpha + beta;
    
   poly   = zeros(dims);
   polyn2 = ones(dims);
   polyn1 = 0.5*(alpha - beta + (alpha + beta + two)*z);
    
   for k = 2:n
    a1 =  two*k*(k + apb)*(two*k + apb - two);
    a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
    a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
    a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);
     
    a2 = a2/a1;
    a3 = a3/a1;
    a4 = a4/a1;

    poly   = (a2 + a3*z).*polyn1 - a4*polyn2;
    polyn2 = polyn1;
    polyn1 = poly;
   end
   jf = poly;
  end

