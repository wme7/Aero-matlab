function f = fvexact(x)

f = exp(-40*(x+.5).^2) + (abs(x-.5)<.25) + eps*rand(size(x));
