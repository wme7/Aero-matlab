function h = numericalFlux(f,df,u,strategy)
% Monotone Riemann fluxes

%% Parameters
nx = length(u); % Grid Size
flx = f(u);     % flux value at every point of the domain
dflx = df(u);   % flux slope at every point of the domain

switch strategy 
    case{1} % Roe Flux
        
        u_ave = (u(1:nx-1) + u(2:nx))/2;    % u @ cells boundaries
        bool_p = df(u_ave) > 0; % if positive
        bool_n = df(u_ave) <= 0;  % if negative
        h = bool_p.*flx(1:nx-1) + bool_n.*flx(2:nx);
        
    case{2} % LF
        
        alpha = max(abs(dflx));
        h = 0.5*(flx(1:nx-1) + flx(2:nx) - alpha*(u(2:nx) - u(1:nx-1)));
        
    case{3} % LLF
        
        fluxmat = [dflx(1:nx-1);dflx(2:nx)];
        beta = max(abs(fluxmat));
        h = 0.5*(flx(1:nx-1) + flx(2:nx) - beta.*(u(2:nx) - u(1:nx-1)));
        
    case{4} % Upwind Flux | Godunov 

        % For dflux is constant along the domain!
        a_p = max(dflx - abs(dflx))/2 == [0]; %#ok<NBRAK> % boolen operator for a>0
        a_n = max(dflx + abs(dflx))/2 == [0]; %#ok<NBRAK> % boolen operator for a<0
        h = a_p*flx(1:nx-1) + a_n*flx(2:nx);

    otherwise 
        error('strategy not suported')
end 