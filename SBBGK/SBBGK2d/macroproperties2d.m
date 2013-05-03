function [r,u1,u2,t,p] = macroproperties2d(n,j_x,j_y,E,nx,ny,theta)
%% Recover Macroscopic Properties
% compute back, fugacity, macroscopic velocities, temperature and pressure.
    % Computing first velocites from the momentum:
    u1 = j_x./n; 
    u2 = j_y./n; 
   
% to compute fugacity, temperature and pressure, we need to rely on the
% distribution fucntion that we where are using: MB, FD, BE.

switch theta
    case{-1} % BE
    % If BE: we apply bisection method to the approx BE distribution Eq.
        r_a = 0.001; r_b = 0.999; tol = 1e-7;
        
        for j = 1:ny
            for i = 1:nx
                psi = @(r_x)2*E(j,i) - BE(r_x,2)/pi*(n(j,i)/BE(r_x,1)).^2 ...
                    - n(j,i).*(u1(j,i).^2+u2(j,i).^2);
                r_p = bisection(psi,r_a,r_b,tol);
                r(j,i) = r_p;
                t(j,i) = n(j,i)/(pi*BE(r_p,1));
                p(j,i) = E(j,i) - 1/2*n(j,i)*(u1(j,i).^2+u2(j,i).^2);
            end
        end
                
    case{1} % FD
    % if FD: we apply bisection method to the approx FD distribution Eq.
        r_a = 0.001; r_b = 0.999; tol = 1e-7;
        
        for j = 1:ny
            for i = 1:nx
                psi = @(r_x)2*E(j,i) - FD(r_x,2)/pi*(n(j,i)/FD(r_x,1)).^2 ...
                    - n(j,i).*(u1(j,i).^2+u2(j,i).^2);
                r_p = bisection(psi,r_a,r_b,tol);
                r(j,i) = r_p;
                t(j,i) = n(i)/(pi*FD(r_p,1));
                p(j,i) = E(j,i) - 1/2*n(j,i)*(u1(j,i).^2+u2(j,i).^2);
            end
        end
    
    case{0} % MB
    % IF MB: the task is much simple.
        t = 2*E./n - (u1.^2 + u2.^2);
        r = n./(pi.*t);
        p = (1/2).*n.*t;
    otherwise 
        error('theta can only be: -1, 0, +1 ');
end