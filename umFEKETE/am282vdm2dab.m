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

 % Build Vandermonde based at nodes
 function vdm = am282vdm2dab(umP, a,b)

 % set the nodes to be equispaced
 umNpts = (umP+1)*(umP+2)/2;
 
 % build the Vandermonde matrix (using the monomials)
 % vdm_ij = (x_i)^j
 dims = size(a);
 vdm = zeros(dims(1)*dims(2),umNpts);

 skip = 1;
 for i=1:umP+1
  for j=1:umP+1 + (1-i)
   vdm(:,skip) = am282jacobi2dab(a,b,i-1,j-1); 
   skip = skip+1;
  end
 end
