
switch method
    case{1} % TVD 0(h^2)
        % Using discrete ordinate method (discrete and constant velocity
        % values in phase-space domain)
        a = v(:,1);        
        % Load initial condition
        f = f0;
%            for tsteps = time
               f_eq = f_equilibrium_1d(r,ux,v,t,theta);
                % initialize variables
                 u_next = zeros(1,nx);
                 u_eq = zeros(1,nx);
                 u = zeros(1,nx);
                 for i = 1:nv
                      % load subcase
                      u_eq(:) = f_eq(i,:);
                      u(:) = f(i,:);
                      % Compute the smoothness factors, r(j), from data, u(j).
                       [r] = theta1d(u,a(i));
                        % Compute the Flux Limiter
                       [phi] = fluxlimiter1d(r,1); % using limiter = 1
                       % Compute TVD Fluxes
                       [F_left,F_right] = TVDflux1d(u,a(i),dtdx,phi);
                       % Compute next time step
                        u_next = u - dtdx*(F_right - F_left) ...
                        + (dt/r_time)*(u_eq-u);
                        % BC
                        u_next(1) = u_next(2);
                        u_next(nx) = u_next(nx-1);
                        % UPDATE info
                         u = u_next;                
                        % Going back to f
                         f(i,:) = u(:);                         
                 end
                   % Compute macroscopic moments
                   [n,j_x,E] = macromoments1d(k,w,f,v);
            
                   % UPDATE macroscopic properties 
                    % (here lies a paralellizing computing chalenge)
                   [r,ux,t,p,yun] = macroproperties1d(n,j_x,E,nx,nv,theta);
%                         [p,yun] = macroproperties1d(n,j_x,E,nx,nv,theta);
%            end
           case{2} % WENO k = 3 i.e. O(h^5) 
    otherwise
          error('Order must be between 1 and 2');  
end

