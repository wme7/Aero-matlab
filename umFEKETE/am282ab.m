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

  function [a,b] = am282ab(x,y)

   dims = size(x);
   Np = dims(1)*dims(2);

   a = zeros(dims);
   b = zeros(dims);
   for n=1:Np  % hack
    if(y(n) ~= 1)
     a(n) = 2*(1+x(n))/(1-y(n))-1;
    else
     a(n) = -1;
    end
   end

   b = y;
