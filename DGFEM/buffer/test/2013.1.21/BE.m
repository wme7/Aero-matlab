function BE = BoseEinstein(z,nu)
%% Bose-Einstein function
% Implementation of the Bose-Einstein Statistical function
% in latex words:
%
% $$F_\nu(z)=\frac{1}{\Gamma(\nu)}\int_{0}^{\infty}\frac{x^{\nu-1}}{z^{-1}e^x-1} 
%   \approx \sum_{l=1}^{\infty}\frac{z^l}{l^\nu}$$
%
% As we can notice this function is of the order $\nu$ 

l   = 1:50; % up to l = 50 to ensure a fair accuaracy
be  = (z.^l) ./ (l.^nu);
BE  = sum(be);
return