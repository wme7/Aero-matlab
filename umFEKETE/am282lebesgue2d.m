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

 % test Lebesgue number

 % set the FINE nodes to be equispaced
 fineP = 100;
 fineNpts = (fineP+1)*(fineP+2)/2;
 fineX = zeros(fineNpts,1);
 fineY = zeros(fineNpts,1);

 skip = 1;
 for i=1:fineP+1
  for j=1:fineP+1 + (1-i)
   fineX(skip) = (-1*(fineP-(i-1))+1*(i-1))/fineP;
   fineY(skip) = (-1*(fineP-(j-1))+1*(j-1))/fineP;
   skip = skip + 1;
  end
 end 
 
 [a,b] = am282ab(X,Y);
 vdmNew  = am282vdm2dab(umP, a,b);

 [fineA,fineB] = am282ab(fineX,fineY);
 vdmFine = am282vdm2dab(umP, fineA,fineB);

 % for each Lagrange interpolating node
 dims = size(vdmNew);
 if(dims(1) ~= dims(2))
  fn = vdmFine*((vdmNew'*vdmNew)\vdmNew');
 else
  fn = vdmFine/vdmNew;
 end

 % plot first node
 fineTRI = delaunay(fineX,fineY);

 trisurf(fineTRI,fineX,fineY,sum(abs(fn)'));

 LebesgueNumber = max(sum(abs(fn')))     




