function h = fluxlimiter(fluxk,fa,fb,a,b)
%% Numerical Flux 'h'
% This are two points Lipschitz continuous monotone flux h(a,b).
% Here, 'fa' and 'fb' are two vector rows in which we have all the
% initial and final value of every element are copied.

%% Compute f'(u) for our data

switch fluxk
    case{1} % Roe with Entropy Fix
        if dfu > 0 
        h = fa;
        elseif dfu < 0
        h = fb;
        else
        h = hllf;
        end
        
    case{2} % (Global) Lax-Friedrichs
        alpha = max(abs(1)); % Is the maximun over the entire region
        h = 1/2*(fa+fb-alpha*(b-a));
        
    case{3} % Local Lax-Friedrichs
        beta  = max(abs(1)); % Is the local maximun in the range min(a,b)<u<max(a,b)
        h = 1/2*(fa+fb-beta*(b-a));
        
    otherwise
        error('Flux limiter case no found, check for available options')
end
