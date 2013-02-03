function phi_x = fluxlimiter1d(r_x,limiter)
%% Compute by case:
switch limiter
    case{1} % Van Leer
        fprintf('\n using Van Leer Limiter \n\n');
        phi_x = (r_x + abs(r_x))./(1 + abs(r_x));
         
    case{2} % Superbee
        fprintf('\n using Superbee Limiter \n\n');
        %phi = max(0,min(2.*r,0),min(r,2));
        phi_x_star = max(min(2.*r_x,0),min(r_x,2));
        phi_x = max(0,phi_x_star);
        
    case{3} % Minmod
        fprintf('\n using Minmod Limiter \n\n');
        phi_x = max(0,min(1,r_x));

    case{4} % koren
        fprintf('\n using Koren Limiter \n\n');
         %phi = max(0,min(2*r,(2/3)*r+1/3,2));
        phi_x_star = min(2*r_x,(2/3)*r_x+1/3);
        phi_x_hat  = min(phi_x_star,2);
        phi_x = max(0,phi_x_hat);

    otherwise
        error('only supported 1,2,3 and 4');
end

return