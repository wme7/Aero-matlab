function int = ModGaussChebyshev(N)
% returns guass-chebyshev integration points and weights transformed to
% work on -inf,inf
  int.n = N;
  int.x = zeros(N,1);
  int.w = zeros(N,1);
  for i=1:N,
    int.x(i) = cos(pi*(i - 0.5)/N);
    int.w(i) = pi/N * 2.0/sqrt(1.0 - int.x(i)^2);
    int.x(i) = log((1.0-int.x(i))/(1.0+int.x(i)));
  end;
end;