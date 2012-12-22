function residue = AdvecResidual(ut,F,dF,S,Ln,Lp,V,D,invM,flux_type)
% Parameters
[np,nx] = size(ut);

%% Residue
% Compute u in element nodes
u = (ut'*V')';

% Compute flux legendre coeficients 'ft' over every element
f = F(u); ft = V\f;

% Compute source legendre coeficients 'st' over every element
s = S(u); st = V\s;

% Contribution of fluxes at cell boundary See Ref. [4]
un = (ut'*Ln)'; % u_{i+1/2}^(-) -> Right u
up = (ut'*Lp)'; % u_{i-1/2}^(+) -> Left u
ub = [up(2:nx);un(1:nx-1)]; 
h  = DGflux1d(F,dF,ub,flux_type); % Evaluate fluxes

% Compute residue function:
residue = zeros(np,nx);
for i = 2:nx-1
    residue(:,i) = ( D'*ft(:,i) - ...           % Volume term
                     h(i)*Ln + h(i-1)*Lp + ...  % flux terms
                     st(:,i) ).*diag(invM);     % Source term
end

% BCs: Periodic BC assuming element 1 and nx are ghost cells.
ut_next(:,1)  = ut(:,nx-1);
ut_next(:,nx) = ut(:,2);