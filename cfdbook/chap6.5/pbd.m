function ye = pbd(t, x, y, side)
%PBD	Prescribes pressure p at outflow boundaries
% Possible values for side: 'lower', 'upper',  'left', 'right'

global geval Re 

ye = 0;
if geval == 1|geval == 7		% Horizontal Poiseuille flow
  if strcmp(side, 'right') == 1		% Outflow
    ye = 0;
  else
    ye = 0;
  end
elseif geval == 2	% Vertical Poiseuille flow
  if strcmp(side, 'upper') == 1		% Outflow
    ye = 0;
  else
    ye = 0;
  end
elseif geval == 3	% Backward facing step
    ye = 0;
elseif geval == 4	% Driven cavity
    ye = 0;
elseif geval == 5|geval == 6	% Uniform flow under angle alpha
  ye = 1; 
else
  error('Wrong value in input for parameter geval')  	
end

