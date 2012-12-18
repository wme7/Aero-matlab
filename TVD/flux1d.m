function h = flux1d(f,df,u,strategy)
% Compute two points Lipschitz continuous monotone fluxes
% based on Cockburn & Shu (1989) Math of Comp. Vol.52, No 186 PP.411-435 
% Coded by Manuel Diaz 2012.12.17

%% Prepare Data
nx = length(u);     % number of data points in our domain
h  = zeros(1,nx-1); % Monotone flux output
flux = f(u);        % flux value at every point of the domain
dflux = df(u);      % flux slope at every point of the domain

%% Compute the flux using 1 the following strategies:
switch strategy
    % assume a = u(i) and b = u(i+1)
    
    case {1} % Godunov
        
        for i = 1:nx-1 % for all the middle points
            % if a <= b choose min f(u) of (a<u<b) 
            % if a >= b choose max f(u) of (a>u>b)
            if u(i) <= u(i+1)
                h(i) = min(flux(i),flux(i+1));
            else % u(i+1) <= u(i)
                h(i) = max(flux(i),flux(i+1));
            end
        end
        
    case {2} % Roe with entropy fix
        
        % Choose f(a) if f'(u) >= 0 for u that belongs [min(a,b),max(a,b)]
        % Choose f(b) if f'(u) <= 0 for u`that belongs [min(a,b),max(a,b)]
        % Otherwise choose LLF flux
        for i = 1:nx-1 % all middle poinst
            u_ave = (u(i) + u(i+1))/2;
            if df(u_ave) > 0 % if positive
                h(i) = f(u(i));
            elseif df(u_ave) <= 0 % if negative
                h(i) = f(u(i+1));
            else
                h(i) = 0;
            end
        end
        
    case {3} % (Global) Lax Friedrichs
        
        % Compute  max(|f'(u)|) 
        % This max value is computed over the entire region of u, that is
        % [inf u(x), sup u(x)] where u(x) is the initial function.
        alpha = max(abs(dflux));
        % Compute Lax-Friedrichs Flux
        for i = 1:nx-1 % for all the middle points
            h(i) = 0.5*( flux(i) + flux(i+1) - alpha *( u(i+1) - u(i) ));
        end
        
    case {4} % Local Lax-Friedrichs
        
        % Similarly for Lax-Friedrichs:
        for i = 1:nx-1 % for all the middle points
            beta = max(abs(dflux(i)),abs(dflux(i+1)));
            h(i) = 0.5*( flux(i) + flux(i+1) - beta *( u(i+1) - u(i) ));
        end
        
    case {5} % Simple Upwind
        
    % For dflux is constant along the domain!
    % if dflux > 0, flux to the left, then h(a,b) = h(a)
    % if dflux > 0, flux to the left, then h(a,b) = h(b)
    % We evaluate professor yang's strategy:
    % % a = (a - |a|)/2
    a = max(dflux - abs(dflux))/2; 
    for i = 1:nx-1 % for all middle points
        if a == 0 % a > 0 
            h(i) = flux(i); % Flux to the left
        else % a < 0 
            h(i) = flux(i+1); % Flux to the right
        end
    end
    
    otherwise
        error('flux strategy not available')
end