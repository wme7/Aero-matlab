function [u,tCells] = limitSolution(u,xgrid,M)
%% Perform WENO Limiting procedure
    w = xgrid.weights';
    nE = xgrid.nElements;
    x = xgrid.nodeCoordinates;
    dx = xgrid.elementSize;
    K = xgrid.solutionDegree;
        
    % Build u_j(x) lagrange polynomials
    ulp = zeros(size(x,2),size(x,1));
    for j = 1:nE
        ulp(j,:) = lagrange(x(:,j),u(:,j)); % interpolation
    end
    
    % Build Cell averages for every E_j
    u_bar = w*u/2;
      
    % Interpolate u and flux values at the boundaries of Ij
    switch xgrid.quadratureType
        case 'LGL'
            u_lbd = u(1,:);
            u_rbd = u(end,:);
        otherwise
            u_lbd = L.lcoef*u;
            u_rbd = L.rcoef*u;
    end
    % Build Numerical fluxes acroos faces
    u_pface = [u_lbd,0]; % + side 
    u_nface = [0,u_rbd]; % - side 

    % Apply Periodic BCs
    %u_nface(1) = u_nface(end); % left BD
    %u_pface(end) = u_pface(1); % right BD
    
    % Apply Neumann BCs
    u_nface(1) = u_pface(1); % left BD
    u_pface(end) = u_nface(end); % right BD
       
    % detect troubled cells
    tcd = TroubleCellDectector(u_bar,'MODminmod',1,u_nface,u_pface,M,dx);
    tCells = tcd.troubledCells;

    for j = tCells
        if ( j~=1 && j~=nE )
            % Build smooth indicators
            for s = 1:K
                dpl = ulp(j-1,:);
                % Derivate 's' times
                for i = 1:s
                    dpl = polyder(dpl);
                end
                integ=dx^(2*s-1)*polyint(conv(dpl,dpl));
                term(s)=polyval(integ,x(1+K,j-1))-polyval(integ,x(1,j-1));
            end
            Beta0 = sum(term);
            
            for s = 1:K
                dpl = ulp(j,:);
                % Derivate 's' times
                for i = 1:s
                    dpl = polyder(dpl);
                end
                integ=dx^(2*s-1)*polyint(conv(dpl,dpl));
                term(s)=polyval(integ,x(1+K,j))-polyval(integ,x(1,j));
            end
            Beta1 = sum(term);
            
            for s = 1:K
                dpl = ulp(j+1,:);
                % Derivate 's' times
                for i = 1:s
                    dpl=polyder(dpl);
                end
                integ=dx^(2*s-1)*polyint(conv(dpl,dpl));
                term(s)=polyval(integ,x(1+K,j+1))-polyval(integ,x(1,j+1));
            end
            Beta2 = sum(term);
            
            % Build Weights
            gamma = [1e-6,0.999998,1e-6]; epsilon = 1e-6;
            w_tilde0 = gamma(1)./(epsilon+Beta0).^2;
            w_tilde1 = gamma(2)./(epsilon+Beta1).^2;
            w_tilde2 = gamma(3)./(epsilon+Beta2).^2;
            wsum = w_tilde0 + w_tilde1 + w_tilde2;
            w0 = w_tilde0./wsum;
            w1 = w_tilde1./wsum;
            w2 = w_tilde2./wsum;
            
            % Modify troubled polynomials
            P0 = ulp(j-1,:);
            P1 = ulp( j ,:);
            P2 = ulp(j+1,:);
            P0_bar = u_bar(j-1);
            P1_bar = u_bar(j);
            P2_bar = u_bar(j+1);
            P0_tilde = P0(K+1) - P0_bar + P1_bar;
            P2_tilde = P2(K+1) - P2_bar + P1_bar;
            u_new = w0*polyval(P0_tilde,x(:,j)) ...
                + w1*polyval(P1,x(:,j)) ...
                + w2*polyval(P2_tilde,x(:,j));
            u(:,j) = u_new; %update info
        end
    end
end    