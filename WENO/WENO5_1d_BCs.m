function u_next = WENO5_1d_BCs(u_next,bc_type,nx)
% Compute boundary conditions for WENO3 in 1D Domain

switch bc_type
    case{1} % Dirichlet BC
        fix_a = 1.0;
        fix_b = 0.0;
        u_next(1)    = fix_a; % left boundary condition (fixed value)
        u_next(2)    = fix_a; % left boundary condition (fixed value)
        u_next(3)    = fix_a; % left boundary condition (fixed value)
        u_next(4)    = fix_a; % left boundary condition (fixed value)
        u_next(nx-3) = fix_b; % right boundary condition (fixed value)
        u_next(nx-2) = fix_b; % right boundary condition (fixed value)
        u_next(nx-1) = fix_b; % right boundary condition (fixed value)
        u_next( nx ) = fix_b; % right boundary condition (fixed value)
        
    case{2} % Neumann BC
        u_next(nx-3) = u_next(nx-4); % right boundary condition
        u_next(nx-2) = u_next(nx-3); % right boundary condition
        u_next(nx-1) = u_next(nx-3); % right boundary condition
        u_next( nx ) = u_next(nx-3); % right boundary condition
        u_next(1)    = u_next(4);    % left boundary condition
        u_next(2)    = u_next(4);    % left boundary condition
        u_next(3)    = u_next(4);    % left boundary condition
        u_next(4)    = u_next(5);    % left boundary condition
        
    case{3} % Periodic BC
        u_next(nx-3) = u_next(3);    % right boundary condition
        u_next(nx-2) = u_next(4);    % right boundary condition
        u_next(nx-1) = u_next(5);    % right boundary condition
        u_next( nx ) = u_next(6);    % right boundary condition
        u_next(1)    = u_next(nx-5); % left boundary condition
        u_next(2)    = u_next(nx-4); % left boundary condition
        u_next(3)    = u_next(nx-3); % left boundary condition        
        u_next(4)    = u_next(nx-2); % left boundary condition        
%    case{4} % Reflective BC
        
    otherwise
        error('case not available')
end