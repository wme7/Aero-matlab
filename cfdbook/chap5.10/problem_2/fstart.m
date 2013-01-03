function ye = fstart(t)		% Initial condition
del = 10^(-15);
if t < 0.5 - del
  y = -1;
elseif t > 0.5 + del
  y = 1;
else
  y = -1 + (y - 0.5 + del)/del;
end
ye = y;
