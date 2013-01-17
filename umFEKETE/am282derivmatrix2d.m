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

 % Build derivative matrix for function held at the nodes
 % [ use orthogonal basis to construct the matrix ]

 function [Dr,Ds] = am282derivmatrix2d(umP,x,y)
 
  dims = size(x);
  umNpts = (umP+1)*(umP+2)/2;
 
  dRhat = zeros(dims(1),umNpts);
  dShat = zeros(dims(1),umNpts);

  % find tensor-product coordinates
  [a,b] = am282ab(x,y); 

  skip = 1;
  for i=1:umP+1
   for j=1:umP+1 + (1-i)
    [dRhat(:,skip),dShat(:,skip)] = am282jacobideriv2dab(a,b,i-1,j-1); 
    skip = skip+1;
   end
  end

  vdm = am282vdm2dab(umP,a,b);
 
  dims = size(vdm);

  if(dims(1) ~= dims(2))
   % f = VDM f~
   % f^ = (VDM'*VDM)^-1 VDM' f
   Dr = dRhat*((vdm'*vdm)\vdm');
   Ds = dShat*((vdm'*vdm)\vdm');
  else
   Dr = dRhat/vdm;
   Ds = dShat/vdm;
  end
