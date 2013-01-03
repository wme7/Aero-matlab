function y = limiter(a,b,limtype)
% Limiter function as defined in Sect. 4.8 and as used in Sect. 10.8
if  b==0,  y=0; else
  if limtype == 2          	% van Albada
    y = (a^2 + a*b)/(a^2 + b^2);
  elseif limtype == 1		% Minmod
    y = max(0, min(a/b,1));
  else				% No limiting, no MUSCL
    y = 0; 
  end   
end
