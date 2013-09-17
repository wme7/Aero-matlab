function res = burgerCPR1d_RHS(u,k,dgl,dgr,coefd,ne,elementdbNodes)
% Compute the Residual:
%
% Globals:  ne, number of elements
%           edbn, elements boundary nodes
% Inputs:   u, data
%           k, polynomial degree,
% Output:   res, residual.

nf = ne+1;
res= zeros(k+1,ne);
qr = zeros(nf);
ql = zeros(nf);

%% Define our Flux function
f = @(w) w.^2/2; %a*w; 
% and the Derivate of the flux function
df = @(w) w; % a*ones(size(w));
    
%% Evaluate flux function in every solution point in domain
q = f(u)

%%  Reimann fluxes for each face
% Mask positive face at x_{i+1/2}^{+}
pFace = elementdbNodes(1,2:end);

% Mask negative face at x_{i+1/2}^{-}
nFace = elementdbNodes(2,1:end-1);

% Compute riemann fluxes (using LF flux spliting)
alfa = max( df(u(pFace)) - df(u(nFace)) );
Nq = 0.5*(alfa*(u(pFace)-u(nFace)) + (q(pFace) + q(nFace)));

%% Fluxes for each face, the left and right state
for i=2:ne
    ql(i)=q(k+1,i-1);
    qr(i)=q(1,i);
end

% Periodic boundary conditions
i=1;
ql(i)=q(k+1,ne);
qr(i)=q(1,1);

i=nf;
ql(i)=q(k+1,ne);
qr(i)=q(1,1);

%% Compute derivative term and correction term
% every element loop
for i=1:ne
    % every solution point loop
    for j=1:k+1
        res(j,i) = 0.;
        % the derivative term
        res(j,i) = res(j,i) + coefd(j,:)*q(:,i);
        % lifting term
        res(j,i) = res(j,i) - alfa(j,1)*(ql(i)-qr(i));
    end    
end