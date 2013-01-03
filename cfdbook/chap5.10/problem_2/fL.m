function ye = fL(t)		% Left boundary value
del = 10^(-15);
if t-fix(t) < 0.5 - del
  y = 1;
elseif t-fix(t) > 0.5 + del
  y = -1;
else
  y = 1 - (y - 0.5 + del)/del;
end
ye = y;
ye = -1;	% Necessary if the velocity is zero
