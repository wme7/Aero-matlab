function ye = exact_solution(t,x,c)

	% Function called: profile
	
yyy = x;
for j = 1:length(x),   yyy(j) = profile(x(j) - c*t); end
ye = yyy;
