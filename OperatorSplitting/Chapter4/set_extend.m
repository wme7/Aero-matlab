function ad=set_extend(Afunc,u,bndcond)
% Fills in ghost cells according to the boundary condition,
% and evaluates Afunc on the extended vector.
%
  S=size(u); N=S(1);
  switch lower(bndcond)
   case{'neumann'}
    uext=[u(1,:);u;u(N,:)];            % "Neumann"
   case{'dirichlet'}
    uext=[zeros(1,S(2));u;zeros(1,S(2))]; % Dirichlet
   case{'periodic'}
    uext=[u(N,:);u;u(1,:)];             % "periodic"
   case {'extrap'}
    uext=[2*u(1,:)-u(2,:);u; 2*u(N,:)-u(N-1,:)];% "extrapolated"
  end;
  ad=feval(Afunc,uext);
