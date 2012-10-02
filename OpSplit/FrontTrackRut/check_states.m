function y=check_states(x,s,t_start,t)
x=x+s.*(t-t_start);
dx=min(diff(x));
if (dx<-1e-6)
	y=0;
	dx
else
	y=1;
end;