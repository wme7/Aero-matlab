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

function [dmodedr, dmodeds] = am282jacobideriv2dab(a,b,ID,JD)

  dims = size(a);

  j2deriv = zeros(dims);

  fa  = am282jacobi1d(a,ID,0,0);
  dfa = am282jacobideriv1d(a,ID,0,0);

  gb  = am282jacobi1d(b,JD,2*ID+1,0);
  dgb = am282jacobideriv1d(b,JD,2*ID+1,0);

  % r-derivative 
  % d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
  dmodedr = dfa.*gb;
  if(ID>0)
   dmodedr = dmodedr.*((0.5*(1-b)).^(ID-1));
  end 
  
  % s-derivative 
  % d/ds = da/ds d/da + db/ds d/db = (2(r+1)/(1-s)^2) d/da + d/db
  %                                = ((1+a)/2)/((1-b)/2) d/da + d/db
  dmodeds = dfa.*(gb.*(0.5*(1+a)));
  if(ID>0)
   dmodeds = dmodeds.*((0.5*(1-b)).^(ID-1)); 
  end 

  tmp = dgb.*((0.5*(1-b)).^ID); 
  if(ID>0)
   tmp = tmp-0.5*ID*gb.*((0.5*(1-b)).^(ID-1));
  end

  dmodeds = dmodeds+fa.*tmp;


   % normalise
  dmodedr = dmodedr*sqrt((ID+0.5)*(ID+JD+1));
  dmodeds = dmodeds*sqrt((ID+0.5)*(ID+JD+1));





