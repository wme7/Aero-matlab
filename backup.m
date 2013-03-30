            % initialize variables
            u_eq = zeros(1,nx);
            u = zeros(1,nx);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                switch f_case
                    case{1} % Relaxation Scheme
                        u_eq(:) = f_eq(i,:);
                        u(:) = f(i,:);
                    case{2} % Euler Limit
                        u_eq(:) = f_eq(i,:);
                        u(:) = f_eq(i,:);
                end
                
                % Compute TVD Fluxes
                [F_left,F_right] = Upwindflux1d(u,a(i,:));

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