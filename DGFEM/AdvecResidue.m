function residue = AdvecResidue(ut,F,dF,S,Ln,Lp,V,D,invM,flux_type)
% Parameters
[np,nx] = size(ut);

%% Residue
% Compute u in element nodes
u = V*ut;

% Compute flux legendre coeficients 'ft' over every element
f = F(u); ft = V\f;

% Compute source legendre coeficients 'st' over every element
s = S(u); st = V\s;
    
% % transform s(x,t) to degress of freedom s(t)_{l,i} for each cell.
% st = zeros(np,nx);
% for l = 0:k             % for all degress of freedom
%     i = l+1;            % Dummy index
%     for j = 1:nx
%         %st(i,j) = (2*l+1)/2.*sum(w(:,j).*s(:,j).*V(:,i));
%         st(i,j) = dx/2*sum(w(:,j).*s(:,j).*V(:,i));
%     end
% end

% Compute fluxes - Inner cells boundaries
up = (ut'*Lp)'; % u_{i-1/2}^(+) -> Left u
un = (ut'*Ln)'; % u_{i+1/2}^(-) -> Right u
uc = [un(1:nx-1);up(2:nx)];
h  = DGflux1d(F,dF,uc,flux_type); % Evaluate fluxes

% Compute fluxes for periodic BC - boundary cells
ub = [un(nx);up(1)];
hb = DGflux1d(F,dF,ub,flux_type);

% Periodic BC @ cell#1
residue(:,1) =  ( D'*ft(:,1) - ...          % Volume term
                h(1)*Ln + hb*Lp  ...        % flux terms
                ).*diag(invM) + st(:,1);    % Source term

% Compute residue function:
for i = 2:nx-1
    residue(:,i) = ( D'*ft(:,i) - ...       % Volume term
                    h(i)*Ln + h(i-1)*Lp ... % flux terms
                ).*diag(invM) + st(:,i);    % Source term
end

% Periodic BC @ cell#nx
residue(:,nx) = ( D'*ft(:,nx) - ...         % Volume term
                hb*Ln + h(nx-1)*Lp  ...     % flux terms
                ).*diag(invM) + st(:,nx);   % Source term
            