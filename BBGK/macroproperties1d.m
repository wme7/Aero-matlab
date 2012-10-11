function [r,u,t,p] = macroproperties1d(n,j_x,E,nx,nv,theta)
%% Recover Macroscopic Properties
% compute back, fugacity, macroscopic velocities, temperature and pressure.
    % Computing first velocites from the momentum:
    u = j_x./n; 
   
% to compute fugacity, temperature and pressure, we need to rely on the
% distribution fucntion that we where are using: MB, FD, BE.

switch theta
    case{-1} % BE
    % If BE: we apply bisection method to the approx BE distribution Eq.
        r_a = 0.0001; r_b = 0.99; 
        for i = 1:nx
        psi = @(x)2*E(i) - BoseEinstein(x,2)/pi*(n(i)/BoseEinstein(x,1)).^2 ...
            - n(i).*(u(i).^2);        
        r_p = bisection(psi,r_a,r_b);
        r(i) = r_p;
        t(i) = n(i)/(pi*BoseEinstein(r_p,1));
        p(i) = E(i) - 1/2*n(i)*(u(i)^2);
        end
        
        
    case{1} % FD
    % if FD: we apply bisection method to the approx FD distribution Eq.
        r_a = 0.0001; r_b = 0.99; 
        for i = 1:nx
        psi = @(x)2*E(i) - FermiDiract(x,2)/pi*(n(i)/FermiDiract(x,1)).^2 ...
            - n(i).*(u(i).^2);        
        r_p = bisection(psi,r_a,r_b);
        r(i) = r_p;
        t(i) = n(i)/(pi*FermiDiract(r_p,1));
        p(i) = E(i) - 1/2*n(i)*(u(i)^2);
        end        
    
    case{0} % MB
    % IF MB: the task is much simple.
        t = 2*E./n - u.^2;
        r = n./(pi.*t);
        p = (1/2).*n.*t;
    otherwise 
        error('theta can only be: -1, 0, +1 ');
end

% Using Discrete Ordinate Method:
    r = repmat(r,nv,1); u = repmat(u,nv,1); 
    t = repmat(t,nv,1); %p = repmat(p,nv,1);