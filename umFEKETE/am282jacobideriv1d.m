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

% ----------------------------------------------------------------
%  jacobd() - derivative of jacobi polynomials - vector version            
%  
%  Get a vector 'poly' of values of the derivative of the n_th order
%  Jacobi polynomial P^(alpha,beta)_N(z) at the np points z.
%
%  To do this we have used the relation 
%
%  d   alpha,beta   1                  alpha+1,beta+1
%  -- P (z)       = -(alpha+beta+n+1) P  (z)
%  dz  n            2                  n-1
%  ----------------------------------------------------------------*/

function jd = jacobd(z, n, alpha, beta)

 jd = zeros(size(z));

  one = 1.0;
  if(n == 0)
    jd(:,:) = 0.0;
  else
    jd = am282jacobi1d(z,n-1,alpha+one,beta+one);
    jd = jd*0.5*(alpha + beta + n + one);
  end

