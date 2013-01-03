function ye = fL(t)
del = 10^(-15);
if t-fix(t) < 0.5 - del
  y = 1;
elseif t-fix(t) > 0.5 + del
  y = -1;
else
  y = 1 - (y - 0.5 + del)/del;
end
%
ye=y;
%ye = -1;     % ye = -1 for pure diffusion problem (no convection) of Fig. 5.17
