function h = DGflux1d(f,df,u,strategy)
% General flux subroutine 
% We use 1 of 4 flux strategies/algorithms optimized for matlab.
% coded by Manuel Diaz, 2012.12.21 (the end of the world)

%% Parameters
nx = length(u); % Grid Size
%x = 1:nx-1;    % Indexes of the flux at the boundary cells (i.e. middle point in my domain)
flx = f(u);     % flux value at every point of the domain
dflx = df(u);   % flux slope at every point of the domain

switch strategy 
    case{1} % Roe Flux
        
        u_ave = (u(1,:) + u(2,:))/2;    % u @ cells boundaries
        bool_p = df(u_ave) > 0; % if positive
        bool_n = df(u_ave) <= 0;  % if negative
        h = bool_p.*flx(1,:) + bool_n.*flx(2,:);
        
    case{2} % LF
        
        alpha = max(max(abs(dflx)));
        h = 0.5*(flx(1,:) + flx(2,:) - alpha*(u(2,:) - u(1,:)));
        
    case{3} % LLF
        
        %fluxmat = [dflx(1:nx-1);dflx(2:nx)];
        beta = max(abs(dflx));
        h = 0.5*(flx(1,:) + flx(2,:) - beta.*(u(2,:) - u(1,:)));
        
    case{4} % Upwind Flux

        % For dflux is constant along the domain!
        a_p = max(max(dflx - abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a>0
        a_n = max(max(dflx + abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a<0
        h = a_p*flx(1,:) + a_n*flx(2,:);

    otherwise 
        error('strategy not suported')
end 