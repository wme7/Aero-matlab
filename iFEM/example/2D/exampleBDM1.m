clc; clear all; close all;
%% mesh information.
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
maxIt = 4;
%% pde information.
BC = '4'; % Change this number to test other boundary conditions.
switch (BC)
case '1' % Homogenous Dirichlet BC.
  pde.f = inline('2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2))','p');
  pde.exactu = inline('sin(pi*p(:,1)).*sin(pi*p(:,2))','p');
  pde.g_D = inline('0','p');pde.g_N = [];bdEdge = [];pde.d=[];
  pde.Du = inline('[pi*cos(pi*p(:,1)).*sin(pi*p(:,2)) pi*sin(pi*p(:,1)).*cos(pi*p(:,2))]','p');
case '2' % Nonhomogenous Dirichlet BC. 
  pde = mixBCdata;pde.d=[];
  bdEdge = setboundary(node,elem,'Dirichlet');
case '3' % Pure Neumann BC.
  pde = mixBCdata;pde.d=[];
  bdEdge = setboundary(node,elem,'Neumann');
case '4' % Mixed BC.
  pde = mixBCdata;pde.d=[];
  bdEdge = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
case '5' % Jump coefficient.
  pde = jumpdata1;
  bdEdge = setboundary(node,elem,'Dirichlet');
  pde.g_N=[];
end
%% Solve 
err = zeros(maxIt,5); N = zeros(maxIt,1);
for i =1:maxIt
    [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
    [u,sigma,M,edge] = PoissonBDM1(node,elem,pde,bdEdge,'direct');
    barycenter = 1/3.*(node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:));
    uexa = pde.exactu(barycenter);
    if (BC=='3') % Pure Neumann BC
    u = u + uexa(end); 
    end
    err(i,1) = getL2error(node,elem,pde.exactu,u);
    err(i,2) = max(abs(u-uexa));
    err(i,3) = getL2errorBDM1(node,elem,pde.Du,sigma,pde.d);
    err(i,4) = getHdiverrorBDM1(node,elem,pde.f,-sigma,[]);
    if (BC~='5')
        sigmaI = faceinterpolate(pde.Du,node,elem,'BDM1'); % have problem for jump case, now way to interpolant.
        err(i,5)=sqrt((sigma-sigmaI)'*M*(sigma-sigmaI));
    end
    N(i) = size(u,1);
end
%% plot
figure(1); hold on; clf;
r1 = showrate(N,err(:,1),2);
hold on
r2 = showrate(N,err(:,2),2,'-*r');
legend('||u-u_h||',['N^{' num2str(r1) '}'],...
       '||u(c_i) - u_h||_{\infty}',['N^{' num2str(r2) '}'],...
       'LOCATION','Best');
hold off
figure(2); hold on; clf;
r3 = showrate(N,err(:,3),2);
hold on
r4 = showrate(N,err(:,4),2,'-*r');
legend('||\sigma - \sigma_h||',['N^{' num2str(r3) '}'],...
    '|\sigma - \sigma_h|_{H(div)}',['N^{' num2str(r4) '}'],...
    'LOCATION','Best');
if (BC~='5')
    r5 = showrate(N,err(:,5),2,'-*r');
    legend('||\sigma - \sigma_h||',['N^{' num2str(r3) '}'],...
        '|\sigma - \sigma_h|_{H(div)}',['N^{' num2str(r4) '}'],...
        '||\sigma_I - \sigma_h||',['N^{' num2str(r5) '}'],...
        'LOCATION','Best');
end