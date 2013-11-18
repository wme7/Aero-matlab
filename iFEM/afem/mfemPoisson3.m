function [err,time,solver,eqn] = mfemPoisson3(node,elem,pde,bdFlag,option,varargin)
%% MFEMPOISSON solve Poisson equation by various mixed finite element methods
%
%   MFEMPOISSON computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% Created by Ming Wang at Nov., 2012.
%
% See also femPoisson Poisson, crack, Lshape 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = cubemesh([-1,1,-1,1,-1,1],1); % default mesh is a square
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = mixBCdata3;                          % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary3(node,elem,'Dirichlet'); 
end

%% Parameters
option = femoption(option);
maxIt = option.maxIt;   maxN = option.maxN; L0 = option.L0;
elemType = option.elemType; refType = option.refType;


%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
end

%% Initialize err
erruL2 = zeros(maxIt,1);   erruIuhL2 = zeros(maxIt,1); 
errsigmaL2 = zeros(maxIt,1); errsigmaHdiv = zeros(maxIt,1); 
errsigmaIsigmah = zeros(maxIt,1); 
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'RT0-P0'  % RT0-P0 mixed FEM 
%             [u,sigma,eqn,info] = Poisson3RT0(node,elem,pde,bdFlag,option);
            [u,sigma,eqn,face] = Poisson3RT0(node,elem,pde,bdFlag,'direct');
            info.solverTime = 0; info.assembleTime = 0; info.itStep = 0;
            info.stopErr = 0; info.flag = 0;
        case 'BDM1-P0' % BDM1-P0 mixed FEM
            [u,sigma,eqn,info] = Poisson3BDM1(node,elem,pde,bdFlag,option);
    end
    % compute error
    tic;
    if isfield(pde,'Du') && isfield(pde,'f')
        if strcmp(elemType,'RT0-P0')
        errsigmaL2(k) = getL2error3RT0(node,elem,pde.Du,sigma);
        errsigmaHdiv(k) = getHdiverror3RT0(node,elem,pde.f,-sigma,[]);
        sigmaI = faceinterpolate3(pde.Du,node,face);
        else
        errsigmaL2(k) = getL2error3BDM1(node,elem,pde.Du,sigma);
        errsigmaHdiv(k) = getHdiverror3BDM1(node,elem,pde.f,-sigma,[]);
        sigmaI = faceinterpolate3(pde.Du,node,elem,'BDM1');
        end
        errsigmaIsigmah(k)=sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
    end
    if isfield(pde,'exactu')
        erruL2(k) = getL2error3(node,elem,pde.exactu,u);
        % interpolation
        uI = Lagrangeinterpolate(pde.exactu,node,elem,'P0');
        area = simplexvolume(node,elem);
        erruIuhL2(k) = sqrt(dot((uI-u).^2,area));
    end
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),erruL2(k), errsigmaL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = size(node,1);
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
       figure(1);  showresult3(node,elem,u);    
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    tic;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
    meshTime(k) = toc;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:k),erruIuhL2(1:k),1,'-*','||u_I-u_h||_{\infty}',...
              N(1:k),erruL2(1:k),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate3(N(1:k),errsigmaIsigmah(1:k),1,'m-+','||\sigma_I - \sigma_h||',...
              N(1:k),errsigmaHdiv(1:k),1,'g-+','|\sigma - \sigma_h|_{div}',... 
              N(1:k),errsigmaL2(1:k),1,'r-*','||\sigma - \sigma_h||');
end

%% Output
err = struct('N',N,'uL2',erruL2(1:k),'uIuhL2',erruIuhL2(1:k),...
             'sigmaL2',errsigmaL2(1:k),'sigmaHdiv',errsigmaHdiv(1:k),...
             'sigmaIsigmahL2',errsigmaIsigmah(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assmble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error
ts = zeros(k,6); ts = char(ts); ts1 = zeros(k,8); ts1 = char(ts1);
display(' #Dof       ||u-u_h||       ||u_I-u_h||    ||sigma-sigma_h||  ||sigma-sigma_h||_{div}');
% format shorte
display([num2str(err.N) ts num2str(err.uL2,'%0.5e') ts num2str(err.uIuhL2,'%0.5e') ...
                ts num2str(err.sigmaL2,'%0.5e')  ts1 num2str(err.sigmaHdiv,'%0.5e')]);
display(' ================================================================================= ');