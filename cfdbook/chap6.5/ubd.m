function ye = ubd(t, j, k, side)
%UBD	Prescribes u at inflow boundaries
% Possible values for side: 'lower', 'upper',  'left', 'right'

global geval u0 J K yu yv yseglen alpha

ye = 0;
if geval == 1		% Horizontal Poiseuille flow to the right
  if strcmp(side, 'left') == 1
    ye = u0(k*(J+1));	% Inflow profile = outflow profile
  else
    ye = 0;
  end
elseif geval == 2	% Vertical Poiseuille flow
    ye = 0;
elseif geval == 3	% Backward facing step
  if strcmp(side, 'left') == 1
    ye = 6*(yu(k) - yseglen(end)).*(yv(end) - yu(k))/(yv(end) - yseglen(end))^3;
  elseif strcmp(side, 'right') == 1
    ye = 6*(yu(k) - yv(1)).*(yv(end) - yu(k))/(yv(end) - yv(1))^3;
  else
    ye = 0;
  end
elseif geval == 4	% Driven cavity
  if strcmp(side, 'upper') == 1
    ye = 1;
  else
    ye = 0;
  end
elseif geval == 5|geval == 6	% Uniform flow under angle alpha
  ye = cos(alpha); 
elseif geval == 7		% Horizontal Poiseuille flow to the left
  if strcmp(side, 'right') == 1
    ye = u0(k*(J+1));	% Inflow profile = outflow profile
  else
    ye = 0;
  end
else
  error('Wrong value in input for parameter geval')  	
end

