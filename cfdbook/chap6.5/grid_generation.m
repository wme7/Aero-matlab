%GRID_GENERATION
global xu yu xv yv dx dy

tic
dx = [];				% Size of primary (i.e. pressure) cell
for s = 1:length(nx)
  dx = [dx, xseglen(s)*ones(1,nx(s))/nx(s)];
end
dy = [];				% Size of primary (i.e. pressure) cell
for s = 1:length(ny)
  dy = [dy, yseglen(s)*ones(1,ny(s))/ny(s)];
end
[DX,DY] = meshgrid(dx,dy);		% Sizes of primary cells
J= length(dx); K = length(dy);

xu = [0,cumsum(dx)];			% x-coordinates of u-nodes
xv = 0.5*(xu(1:end-1)+xu(2:end));	% x-coordinates of v-nodes
yv = [0,cumsum(dy)];			% y-coordinates of v-nodes
yu = 0.5*(yv(1:end-1)+yv(2:end));	% y-coordinates of u-nodes
[XU,YU] = meshgrid(xu,yu);
[XV,YV] = meshgrid(xv,yv);
[XP,YP] = meshgrid(xv,yu);		% Coordinates of pressure nodes

DXU = XU;			% Size of finite volumes for u; preallocation.
DXU(:,1) = DX(:,1)/2;  		DXU(:,2:J) = (DX(:,1:J-1) + DX(:,2:J))/2;
DXU(:,J+1) = DX(:,J)/2;  	DYU = [DY,dy'];			
DYV = XV;			% Size of finite volumes for v; preallocation.
DYV(1,:) = DY(1,:)/2;   	DYV(2:K,:) = (DY(1:K-1,:) + DY(2:K,:))/2;
DYV(K+1,:) = DY(K,:)/2;		DXV = [DX;dx];
volu = DXU.*DYU; vol = 1./volu';  n = (J+1)*K;
invvolu = spdiags(vol(:), 0, n, n);		% (u-volume)^(-1) 
volv = DXV.*DYV; vol = 1./volv';  n = J*(K+1);
invvolv = spdiags(vol(:), 0, n, n);		% (v-volume)^(-1)

jseg = [0,cumsum(nx)];	% j-indices of horizontal segment boundaries
kseg = [0,cumsum(ny)];	% k-indices of vertical segment boundaries
hx = min(dx);	% Required for stability estimate for  
hy = min(dy);	% Adams-Bashforth-Crank-Nicolson scheme

clear vol 
tijd = toc; disp(['grid_generation time = ',num2str(tijd)])
