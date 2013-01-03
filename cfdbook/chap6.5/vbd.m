function ye = vbd(t, j, k, side)
%VBD	Prescribes v at inflow boundaries
% Possible values for side: 'lower', 'upper',  'left', 'right'

global geval v0 J K alpha

ye = 0;
if geval == 1|geval == 7		% Horizontal Poiseuille flow 
  ye = 0;
elseif geval == 2	% Vertical Poiseuille flow
  if strcmp(side, 'lower') == 1
    ye = v0(j+K*J);	% Inflow profile = outflow profile
  else
    ye = 0;
  end
elseif geval == 3	% Backward facing step
  ye = 0;
elseif geval == 4	% Driven cavity
  ye = 0;
elseif geval == 5|geval == 6	% Uniform flow under angle alpha
  ye = sin(alpha); 
else
  error('Wrong value in input for parameter geval')  	
end

