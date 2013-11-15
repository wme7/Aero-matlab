function y = f(p34)	% Basic shock tube relation equation (10.51)
global PRL  CRL MACHLEFT  gamma

wortel = sqrt(2*gamma*(gamma-1+(gamma+1)*p34));
yy = (gamma-1)*CRL*(p34-1)/wortel;
yy = (1 + MACHLEFT*(gamma-1)/2-yy)^(2*gamma/(gamma-1));
y = yy/p34 - PRL;
