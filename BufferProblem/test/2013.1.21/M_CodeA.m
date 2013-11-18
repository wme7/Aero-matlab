function [rho,u,p,fl,Ml_eq] = ModSBBGK()
%function [rho_our p_out t_out]=M_codeA(rho_in,p_in,t_in)
 %switch method
 %   case{1} % TVD 0(h^2)
        % Using discrete ordinate method (discrete and constant velocity
        % values in phase-space domain)
        a = v(:,1);        
        % Load initial condition
%         if()
        f = f0;
%         end
% %            for tsteps = time
               f_eq = f_equilibrium_1d(rho_0,ux,v,t,RRR);
                % initialize variables                
                 u_next = zeros(1,nx);
                 u_eq = zeros(1,nx);
                 u_l=zeros(1,nx);
                 u = zeros(1,nx);
                 for i = 1:nv
                      % load subcase
                      u_eq(:) = f_eq(i,:);
                      u(:) = f(i,:).*h;
                      u_l(:)=(1.-h).*f(i,:);
                      % Compute the smoothness factors, r(j), from data, u(j).
                       [r] = theta1d(u,a(i));
                        % Compute the Flux Limiter
                       [phi] = fluxlimiter1d(r,1); % using limiter = 1
                       % Compute TVD Fluxes
                       [F_left,F_right] = TVDflux1d(u,a(i),dtdx,phi);
                        [FL_left,FL_right] = TVDflux1d(u_l,a(i),dtdx,phi);
                       % Compute next time step
                        u_next= u - dtdx*(F_right - F_left).*h -dtdx*(FL_right - FL_left).*h + (dt/r_time).*h.*(u_eq-f(i,:));
                            
                     
                        % BC.4
                        u_next(1) = u_next(2);
                        u_next(nx) = u_next(nx-1);
                        % UPDATE info
                         u = u_next;       
                         
%                          for i=1:nx
%                              u
%                          end
                        % Going back to f
                         f(i,:) = u(:);                         
                 end
                   % Compute macroscopic moments
                   
                   [n,j_x,E] = macromoments1d(k,w,f,v);
                 
                   % UPDATE macroscopic properties 
                     [ux,t,p,yun] = macroproperties1d(n,j_x,E,nx,nv,theta);
%                      [r,ux,t,p,yun] = macroproperties1d(n,j_x,E,nx,nv,theta);
%                         [p,yun] = macroproperties1d(n,j_x,E,nx,nv,theta);
%            end
%           case{2} % WENO k = 3 i.e. O(h^5) 
%   otherwise
%          error('Order must be between 1 and 2');  
% end

