function res = burgers1d_RHS(u,k,dgl,dgr,coefd,ne,dbNodes)
% Compute the Residual:
%
% Globals:  ne, number of elements
%           edbn, elements boundary nodes
% Inputs:   u, data
%           k, polynomial degree,
% Output:   res, residual.

nf = ne+1;
res= zeros(k+1,ne);
ur = zeros(1,nf);
ul = zeros(1,nf);

%% Define our Flux function
f = @(w) w.^2/2; %a*w; 
% and the Derivate of the flux function
df = @(w) w; % a*ones(size(w));
    
%% Evaluate flux & derivative flux function in every solution point in domain
q = f(u);

%%  Reimann fluxes for each face
% Mask positive faces (x_{i+1/2}^{+}) & negative face (x_{i+1/2}^{-})
%pFace = dbNodes(1,2:end);
%nFace = dbNodes(2,1:end-1);

% Compute Riemann fluxes (using LF flux spliting)
%alfa = max( df(u(pFace)) - df(u(nFace)) );
%Nq = 0.5*(alfa*(u(pFace)-u(nFace)) + (q(pFace) + q(nFace)));

%% Fluxes for each face, the left and right state
for i=2:ne
    ul(i)=u(k+1,i-1);
    ur(i)=u(1,i);
end

% Periodic BCs
i=1;
ul(i)=u(k+1,ne);
ur(i)=u(1,1);

i=nf;
ul(i)=u(k+1,ne);
ur(i)=u(1,1);

ql = f(ul); qr = f(ur);

% Compute Riemann fluxes (using LF flux spliting)
alfa = max( df(ur) - df(ul) );
Nq = 0.5*(alfa*(ur - ul) + (qr + ql));

%% Compute derivative term and correction term
% Every element loop
for i=1:ne
    % Every solution point loop
    for j=1:k+1
        res(j,i) = 0.;
        % the derivative term
        res(j,i) = res(j,i) + coefd(j,:)*u(:,i);
        % the correction terms
        res(j,i) = res(j,i) + (Nq(i)-ql(i))*dgl(j,1) + (Nq(i)-qr(i))*dgr(j,1);
    end    
end