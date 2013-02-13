function [vp,vn] = WENO_fluxsplit(u,f,df,strategy)
% WENO flux spliting subroutine.
% OUTPUT:
%   * vp: positive flux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative flux  v^{-}, which corresponds to f_{i+1/2}^{+}

switch strategy
    case{1} % Godunov (non-conservative)
        v = f(u);
        vp = 0.5*(v + abs(v)); %flux^{+}
        vn = 0.5*(v - abs(v)); %flux^{-}
    case{2} % Local Lax-Friedrichs
        v = f(u); alpha = abs(df(u));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    case{3} % (Global) Lax-Friedrichs
        v = f(u); alpha = max(abs(df(u)));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    otherwise
        error('only cases 1 and 2 are available')
end