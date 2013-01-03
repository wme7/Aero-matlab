function y = hinit(x,geval)
if geval == 1	% Piecewise linear initial depth
  del = 0.003;
  if x < 0.5 - del,    				y1 = 1;
  elseif (x >= 0.5 - del & x <= 0.5 + del),   	y1 = 1- (x-0.5+del)/(2*del);
  else,    					y1 = 0;
  end
end
if geval == 2	% Sinusoidal initial depth
  y1 = sin(pi*x);
end
if geval == 3 	% Hyperbolic tangent initial depth 
  del=0.05;  y1 = tanh((x-0.5)/del);
end
y=y1;
