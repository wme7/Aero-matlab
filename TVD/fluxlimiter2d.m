function [phi_x,phi_y] = fluxlimiter2d(r_x,r_y,limiter)
%% Compute by case:
switch limiter
    case{1} % Van Leer
        %fprintf('\n using Van Leer Limiter \n\n');
        phi_x = (r_x + abs(r_x))./(1 + abs(r_x));
        phi_y = (r_y + abs(r_y))./(1 + abs(r_y));
        
    case{2} % Superbee
        %fprintf('\n using Superbee Limiter \n\n');
        %phi = max(0,min(2.*r,0),min(r,2));
        phi_x_star = max(min(2.*r_x,0),min(r_x,2));
        phi_x = max(0,phi_x_star);
        phi_y_star = max(min(2.*r_y,0),min(r_y,2));
        phi_y = max(0,phi_y_star);
        
    case{3} % Minmod
        %fprintf('\n using Minmod Limiter \n\n');
        phi_x = max(0,min(1,r_x));
        phi_y = max(0,min(1,r_y));
        
    case{4} % koren
        %fprintf('\n using Koren Limiter \n\n');
         %phi = max(0,min(2*r,(2/3)*r+1/3,2));
        phi_x_star = min(2*r_x,(2/3)*r_x+1/3);
        phi_x_hat  = min(phi_x_star,2);
        phi_x = max(0,phi_x_hat);
        phi_y_star = min(2*r_y,(2/3)*r_y+1/3);
        phi_y_hat  = min(phi_y_star,2);
        phi_y = max(0,phi_y_hat);
        
    otherwise
        error('only supported 1,2,3 and 4');
end

return