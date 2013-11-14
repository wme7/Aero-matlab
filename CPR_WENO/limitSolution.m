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
            
    % detect troubled cells
    tcd = TroubleCellDectector(u_bar,'MODminmod',1,u(1,:),u(1+K,:),M,dx);
    tCells = tcd.troubledCells;

    term0 = zeros(1,size(ulp,2));
    term1 = zeros(1,size(ulp,2));
    term2 = zeros(1,size(ulp,2));
    for j = tCells
        if ( j~=1 && j~=nE )
            % Build smooth indicators
            for s = 1:K
                dpl0 = ulp(j-1,:);
                dpl1 = ulp(j,:);
                dpl2 = ulp(j+1,:);
                % Derivate 's' times
                for i = 1:s
                    dpl0 = polyder(dpl0);
                    dpl1 = polyder(dpl1);
                    dpl2 = polyder(dpl2);
                end
                % Do: dx^(2s-1)*integrate(dp^2,a,b)
                integ0 = dx^(2*s-1)*polyint(conv(dpl0,dpl0));
                term0(s) = polyval(integ0,x(1+K,j-1))-polyval(integ0,x(1,j-1));
                integ1 = dx^(2*s-1)*polyint(conv(dpl1,dpl1));
                term1(s) = polyval(integ1,x(1+K,j))-polyval(integ1,x(1,j));
                integ2 = dx^(2*s-1)*polyint(conv(dpl2,dpl2));
                term2(s) = polyval(integ2,x(1+K,j+1))-polyval(integ2,x(1,j+1));
            end
            Beta0 = sum(term0);
            Beta1 = sum(term1);
            Beta2 = sum(term2);
            
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
            P0_bar = zeros(size(P0));
            P1_bar = zeros(size(P1));
            P2_bar = zeros(size(P2));
            P0_bar(K+1) = u_bar(j-1);
            P1_bar(K+1) = u_bar(j);
            P2_bar(K+1) = u_bar(j+1);
            P0_tilde = P0 - P0_bar + P1_bar;
            P2_tilde = P2 - P2_bar + P1_bar;
            u_new = w0*P0_tilde + w1*P1 + w2*P2_tilde;
            u(:,j) = polyval(u_new,x(:,j)); %update info
        end
    end
    
end    