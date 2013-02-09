function coef = c_rj(r,j,k)
% Generate ENO/WENO stencils coeficients
%**************************************************************************
% Based on:
% Chi-Wang Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% Coded by Manuel Diaz, 02.09.2012, NTU Taiwan.
% Compare with Eqs. 2.20 & 2.21 for uniform grids 
%**************************************************************************

%% Compute coeficient c_rj
sum2 = 0;
for m = j+1:k  
    %% Numerator
    sum = 0; % dummy variable to acumulate the sumation result.
    for l = 0:k
        if l ~= m
            % q sequence product
            prod = 1; % dummy variable to acumulate the product result.
            for q = 0:k
                if q ~= l && q ~= m
                    prod = prod*(r-q+1);
                end
            end
            sum = sum + prod;
        end
    end
    
    %% Denominator
    prod2 = 1; % dummy variable to acumulate the product result.
    for l = 0:k
        if l ~= m
            prod2 = prod2*(m-l);
        end
    end
    %% Total
    sum2 = sum2 + sum/prod2;
end
coef = sum2;