function FD = FermiDiract(z,nu)
%% Fermi-Diract function
% Implementation of the Fermi-Diract Statistical function
% in latex words:
%
% $$F_\nu(z)=\frac{1}{\Gamma(\nu)} \int_{0}^{\infty}\frac{x^{\nu-1}}{z^{-1}e^x+1}
%   \approx \sum_{l=1}^{\infty}(-1)^{l-1}\frac{z^l}{l^\nu}$$
%
% As we can notice this function is of the order $\nu$ 

l   = 1:50; % up to l = 50 to ensure a fair accuaracy
fd  = (-1).^(l-1) .* (z.^l) ./ (l.^nu);
FD  = sum(fd);
return