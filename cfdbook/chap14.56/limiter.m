function y = limiter(a,b,limtype)	% Called by extrap
% Theory in Section 4.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

if  abs(b) < 1/10^9,  y=0;
else
  if limtype == 3	      	% superbee
    r = a/b;    	y = max(0, max(min(2*r,1), min(r,2)));         
  elseif limtype == 2         	% van Albada
    y = (a^2 + a*b)/(a^2 + b^2);
  elseif limtype == 1		% Minmod (Roe(1986))
    y = max(0, min(a/b,1));
  elseif limtype == 4		% PL limiter
    r = a/b;
    if r<0,      	y = 0;
    elseif r<0.4,      	y = 2*r;
    elseif r<4,		y = (2+r)/3;
    else,		y = 2;
    end
  elseif limtype == 5	% Chakravarty/Osher limiter (Kr\"oner (1997) p. 109)
    r = a/b;		y = max(0, min(r,1.1));
  else			       % no MUSCL: first order upwind
    y = 0; 
  end   
end
