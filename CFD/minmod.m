
function mm = minmod(a,b)

     mm = zeros(size(a));
     mm =    ((abs(a)<=abs(b)).*(a.*b>0)).*a;
     mm = mm+((abs(b)<abs(a)).*(a.*b>0)).*b;
