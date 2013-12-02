function table = ENO_stencils()
%% Generate Table of ENO/WENO Stencils coeficients
%**************************************************************************
% Based on:
% Chi-Wang Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.09.2012, NTU Taiwan.
% Compare with Eqs. 2.20 & 2.21 for uniform grids!
%**************************************************************************

% Build Table 2.21.
    c = zeros(8,7,7);
    for k = 1:7
        for i = 1:k+1 % dummy index 
            r = i-2; % for range -1 to k-1
            for n = 1:7 % dummy index 
                j = n-1; % for range 0 to 6
                c(i,n,k) = c_rj(r,j,k);
            end
        end
    end
    table = c;
end


function coef = c_rj(r,j,k)
% Generate ENO/WENO stencils coeficients
%**********************************************************************
% Based on:
% Chi-Wang Shu's Lectures notes, 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% Coded by Manuel Diaz, 02.09.2012, NTU Taiwan.
% Compare with Eqs. 2.20 & 2.21 for uniform grids 
%**********************************************************************

% Compute coeficient c_rj
    sum2 = 0;
    for m = j+1:k  
        % Numerator
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

        % Denominator
        prod2 = 1; % dummy variable to acumulate the product result.
        for l = 0:k
            if l ~= m
                prod2 = prod2*(m-l);
            end
        end
        % Total
        sum2 = sum2 + sum/prod2;
    end
    coef = sum2;
end

