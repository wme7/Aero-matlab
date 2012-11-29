function [uun,uup] = WENO3_1D_flux(u)
% Constants
d0 = 1/10; d1 = 6/10; d2 = 3/10;
epsilon = 1E-7;

%% Compute WENO3 1d Fluxes:
% Center value
i = 3;
un = u; up = u;

%% Positive Flux
% Nonlinear Smoothness Indicators
IS0p = 13/12*(up(i-2)-2*up(i-1)+up(i))^2 + 1/4*(up(i-2)-4*up(i-1)+3*up(i))^2;
IS1p = 13/12*(up(i-1)-2*up(i)+up(i+1))^2 + 1/4*(up(i-1)-up(i+1))^2;
IS2p = 13/12*(up(i)-2*up(i+1)+up(i+2))^2 + 1/4*(3*up(i)-4*up(i+1)+up(i+2))^2;

% Negative Flux
% Nonlinear Smoothness Indicators
IS0n = 13/12*(un(i+1)-2*un(i+2)+un(i+3))^2 + 1/4*(3*un(i+1)-4*un(i+2)+un(i+3))^2;
IS1n = 13/12*(un(i)-2*un(i+1)+un(i+2))^2 + 1/4*(un(i)-un(i+2))^2;
IS2n = 13/12*(un(i-1)-2*un(i)+un(i+1))^2 + 1/4*(un(i-1)-4*un(i)+3*un(i+1))^2;
%end

%% Positive Flux
% Non-normilized Stencil weights
alpha0p = d0./(epsilon + IS0p).^2;
alpha1p = d1./(epsilon + IS1p).^2;
alpha2p = d2./(epsilon + IS2p).^2;

% Weigths
w0p = alpha0p ./ (alpha0p + alpha1p + alpha2p);
w1p = alpha1p ./ (alpha0p + alpha1p + alpha2p);
w2p = alpha2p ./ (alpha0p + alpha1p + alpha2p);

% Negative Flux
% Non-normilized Stencil weights
alpha0n = d0./(epsilon + IS0n).^2;
alpha1n = d1./(epsilon + IS1n).^2;
alpha2n = d2./(epsilon + IS2n).^2;

% Weigths
w0n = alpha0n ./ (alpha0n + alpha1n + alpha2n);
w1n = alpha1n ./ (alpha0n + alpha1n + alpha2n);
w2n = alpha2n ./ (alpha0n + alpha1n + alpha2n);

%for i = 3:nx-3
% Positive Flux
% Numerical Flux: u_i+1/2 (+)
uup(i) = w0p*(2/6*up(i-2)-7/6*up(i-1)+11/6*up(i))+...
    w1p*(-1/6*up(i-1)+5/6*up(i)+2/6*up(i+1))+...
    w2p*(2/6*up(i)+5/6*up(i+1)-1/6*up(i+2));

% Negative Flux
% Numerical Flux: u_i+1/2 (-)
uun(i) = w0n*(-1/6*un(i-1)+5/6*un(i)+2/6*un(i+1))+...
    w1n*(2/6*un(i)+5/6*un(i+1)-1/6*un(i+2))+...
    w2n*(11/6*un(i+1)-7/6*un(i+2)+2/6*un(i+3));
end