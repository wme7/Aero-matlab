function flx = flux(u,du,strategy)
% Compute the flux using 1 the 

switch strategy
    case {1} % Godunov
        
    case {2} % Roe with entropy fix
        
    case {3} % (Local) Lax Friedrichs
        
    case {4} % Global Lax Friedrichs
        
    otherwise
        error('flux strategy not available')
end