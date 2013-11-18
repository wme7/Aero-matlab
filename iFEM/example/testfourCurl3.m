%% Test fourCurl3
%
% Lin Zhong, June 2013.



close all; clear all;clc;
pde = fourCurl3data;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet');

maxIt = 3;
errwIwhL2 = zeros(maxIt,1); errwIwhHcurl = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);errwwhHcurl = zeros(maxIt,1);
erruIuhL2 = zeros(maxIt,1); erruIuhHcurl = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);erruuhHcurl = zeros(maxIt,1);
assembleTime = zeros(maxIt,1);
N = zeros(maxIt,1);
option = [];

for i = 1:maxIt
    [w,u,eqn,inf] = fourCurl3(node,elem,bdFlag,pde,option);
    assembleTime(i) = inf.assembleTime;
    uI = edgeinterpolate(pde.exactu,node,eqn.edge);
    wI = edgeinterpolate(pde.curlcurlu,node,eqn.edge);
    errwL2(i) = getL2error3NE(node,elem,pde.curlcurlu,w);
     erruL2(i) = getL2error3NE(node,elem,pde.exactu,u);   
     errwIwhL2(i) = sqrt((w-wI)'*eqn.Me*(w-wI));
     erruIuhL2(i) = sqrt((u-uI)'*eqn.Me*(u-uI)); 
     errwwhHcurl(i) = getHcurlerror3NE(node,elem,pde.curlcurlcurlu,w);
     erruuhHcurl(i) = getHcurlerror3NE(node,elem,pde.curlu,u);   
      errwIwhHcurl(i) = sqrt((w-wI)'*eqn.B*(w-wI));
     erruIuhHcurl(i) = sqrt((u-uI)'*eqn.B*(u-uI));     
     N(i) = size(node,1);
     N(i) = N(i)^(1/3);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end

h = 1./N(1:maxIt);
subplot(2,2,1)
showrateh2(h, errwIwhL2(1:maxIt), 1 , 'k-+', '||w_I-w_h||',...
                   h, erruIuhL2(1:maxIt), 1 , 'r-*', '||u_I-u_h||' );
subplot(2,2,2)
showrateh2(h, errwL2(1:maxIt), 1, 'k-+', '||w-w_h||',...
                   h, erruL2(1:maxIt), 1 , 'r-*', '||u-u_h||' );                                     
subplot(2,2,3)               
showrateh2(h, errwIwhHcurl(1:maxIt), 1 , 'b-+', '||w_I-w_h||_1',...
                   h, erruIuhHcurl(1:maxIt), 1 , 'g-*', '||u_I-u_h||_1' );
subplot(2,2,4)               
showrateh2(h, errwwhHcurl(1:maxIt), 1, 'b-+', '||w-w_h||_1',...
                   h, erruuhHcurl(1:maxIt), 1 , 'g-*', '||u-u_h||_1' );
                                                            
               
% Output
err = struct('N', N,...
                    'wL2',errwL2(1:maxIt),'uL2',erruL2(1:maxIt),...
                    'wIwhL2',errwIwhL2(1:maxIt), 'uIuhL2',erruIuhL2(1:maxIt),...
                    'wIwhHcurl',errwIwhHcurl(1:maxIt),'uIuhHcurl',erruIuhHcurl(1:maxIt),...
                    'wwhHcurl',errwwhHcurl(1:maxIt),'uuhHcurl',erruuhHcurl(1:maxIt));
time = struct('N',N,'assmble',assembleTime(1:maxIt));
       
%% Display error on screen
ts = zeros(maxIt,3); ts = char(ts);
% error of u
fprintf('===========================================================\n');
display('   N  ||w_I-w_h||      ||w-w_h||       ||w_I-w_h||_1        ||w-w_h||_1');  
display([num2str(err.N)  ts num2str(err.wIwhL2,'%0.5e') ts num2str(err.wL2,'%0.5e') ts num2str(err.wIwhHcurl,'%0.5e')  ts num2str(err.wwhHcurl,'%0.5e')]);

fprintf('===========================================================\n');
display('   N    ||u_I-u_h||      ||u-u_h||       ||u_I-u_h||_1        ||u-u_h||_1     ');
display([num2str(err.N) ts num2str(err.uIuhL2,'%0.5e')  ts num2str(err.uL2,'%0.5e')  ts num2str(err.uIuhHcurl,'%0.5e') ts num2str(err.uuhHcurl,'%0.5e')]);
              