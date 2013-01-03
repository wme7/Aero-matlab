function ye = exact_solution(t, x, y, side)
%EXACT_SOLUTION	
% Possible values for side: 'lower', 'upper',  'left', 'right'

global geval xu yu xv yv 

ye = 0;
if geval == 1		% Horizontal Poiseuille flow
  if strcmp(side, 'right') == 1
    ye = 6*(y - yv(1)).*(yv(end) - y)/(yv(end) - yv(1))^3;
  else
    ye = 0;
  end
elseif geval == 2	% Vertical Poiseuille flow
  if strcmp(side, 'upper') == 1
    ye = 6*(x - xu(1)).*(xu(end) - x)/(xu(end) - xu(1))^3;
  else
    ye = 0;
  end
elseif geval == 7		% Horizontal Poiseuille flow
  if strcmp(side, 'right') == 1
    ye = -6*(y - yv(1)).*(yv(end) - y)/(yv(end) - yv(1))^3;
  else
    ye = 0;
  end
else
  error('No exact solution for this case')  	
end

