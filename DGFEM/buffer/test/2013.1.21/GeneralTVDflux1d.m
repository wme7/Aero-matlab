%function [flux] = GeneralTVDflux1d(f,df,u,dtdx,phi_x,strategy)
function [h] = GeneralTVDflux1d(f,df,u,dtdx,phi_x,strategy)
% General TVD flux routine 
% We use 1 of 4 flux strategies/algorithms optimized for matlab.

%% Parameters
nx = length(u); % Grid Size
x = 1:nx-1;     % Indexes of the flux at the boundary cells (i.e. middle point in my domain)
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
        
    case{4} % Upwind Flux

        % For dflux is constant along the domain!
        a_p = max(dflx - abs(dflx))/2 == [0]; % boolen operator for a>0
        a_n = max(dflx + abs(dflx))/2 == [0]; % boolen operator for a<0
        h = a_p*flx(1:nx-1) + a_n*flx(2:nx);

    otherwise 
        error('strategy not suported')
end 

%% Compute TVD fluxes
% Note:
% s = Left flux & r = Right flux
% l = low flux  & h = high flux

% Initalize Arrays
% flux = zeros(1,nx); f_l = zeros(1,nx); f_h = zeros(1,nx);
%     
% for j = x % all middle points
%     % Compute fluxes for TVD
%     f_l(j)  = h(j); % Upwind flux
%     f_h(j)  = 0.5*(flx(j)+flx(j+1))-0.5*dflx(j)*dtdx*(flx(j+1)-flx(j));
%     flux(j) = f_l(j) + phi_x(j)*( f_h(j) - f_l(j) );
% end

return