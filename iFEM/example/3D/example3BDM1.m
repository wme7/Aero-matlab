clc;clear;
%profile on
%% mesh information
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];  % nodes
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7]; % elements
elem = label3(node,elem);
%[node,elem] = bisect3(node,elem,'all');
%[node,elem] = bisect3(node,elem,'all');
%[node,elem] = uniformbisect3(node,elem);
% [node,elem,~,bdFace] = uniformbisect3(node,elem,[],bdFace);
maxIt = 5;
%% pde information.
BC = '1'; % Change this number to test other boundary conditions.
switch (BC)
case '1' % Homogenous Dirichlet BC.
  pde.f = inline('3*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2)).*sin(pi*p(:,3))','p');
  pde.exactu = inline('sin(pi*p(:,1)).*sin(pi*p(:,2)).*sin(pi*p(:,3))','p');
  pde.Du = inline('[pi*cos(pi*p(:,1)).*sin(pi*p(:,2)).*sin(pi*p(:,3)),pi*sin(pi*p(:,1)).*cos(pi*p(:,2)).*sin(pi*p(:,3)),pi*sin(pi*p(:,1)).*sin(pi*p(:,2)).*cos(pi*p(:,3))]','p');
  pde.g_D = [];pde.g_N = [];bdFace = [];
case '2' % Nonhomogenous Dirichlet BC. 
  pde = mixBCdata3;
  bdFace = setboundary3(node,elem,'Dirichlet');
case '3' % Pure Neumann BC.
  pde = mixBCdata3;
  bdFace = setboundary3(node,elem,'Neumann');
case '4' % Mixed BC.
  pde = mixBCdata3;
  bdFace = setboundary3(node,elem,'Dirichlet','all','Neumann','y==-1');
end
%% Solve and plot
err = zeros(maxIt,4); N = zeros(maxIt,1);
for i =1:maxIt
    [node,elem,~,bdFace] = uniformbisect3(node,elem,[],bdFace);
%    [node,elem] = bisect3(node,elem,'all');
    barycenter = 1/4.*(node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:)+node(elem(:,4),:));
    uexa = pde.exactu(barycenter);
    [u,sigma] = Poisson3BDM1(node,elem,pde,bdFace);
    if(BC == '3') % Pure Neumann BC
    u = u + uexa(end);
    end
    err(i,1) = getL2error3(node,elem,pde.exactu,u,2);
    err(i,2) = max(abs(u-uexa));
  %  err(i,3) = getL2Err3BDM1(node,elem,pde.Du,sigma,[]);
  %  err(i,4) = getHdivErr3BDM1(node,elem,pde.f,-sigma,[]);
    N(i) = length(u) + length(sigma) ;
end

save('plotdata.mat','N','err');

figure(1); hold on; clf;
r1 = showrate(N,err(:,1),2);
hold on
r2 = showrate(N,err(:,2),2,'-*r');
legend('||u-u_h||',['N^{' num2str(r1) '}'],...
      '||u(c_i) - u_h||_{\infty}',['N^{' num2str(r2) '}'],...
      'LOCATION','Best')
hold off;
% figure(2); hold on; clf;
% r3 = showrate(N,err(:,3),2);
% hold on
% r4 = showrate(N,err(:,4),2,'-*r');
% legend('||\sigma-\sigma_h||',['N^{' num2str(r3) '}'],...
%       '||\sigma-\sigma_h||_{H(div)}',['N^{' num2str(r4) '}'],...
%       'LOCATION','Best')
% hold off;
%profile viewer
