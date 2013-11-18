function [err,time,solver,eqn] = femStokesHdiv(node,elem,pde,bdFlag,option,varargin)
%% FEMPOISSON solve Poisson equation by various finite element methods
% Created by Ming Wang, at Nov., 2012.
%

%% Check input arguments
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = squaremesh([0,1,0,1],0.25);  % default mesh is a square
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
   pde = StokesZulehnerdata;                   % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end
% default elemType
if ~isfield(option,'elemType')
    option.elemType = 'RT0-P0';
end

%% Parameters
option = femoption(option);
elemType = option.elemType;
maxIt = option.maxIt;
maxN = option.maxN;
L0 = option.L0;
refType = option.refType;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
erruL2 = zeros(maxIt,1); erruIuhH1 = zeros(maxIt,1); erruInf = zeros(maxIt,1);
errpL2 = zeros(maxIt,1); errpIphL2 = zeros(maxIt,1); errpInf = zeros(maxIt,1);
errpL2re = zeros(maxIt,1);
errwL2 = zeros(maxIt,1); errwIwh = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'RT0-P0'   % RT0-P0 mixed FEM 
            [u,p,w,edge,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option);
        case 'BDM1B-P0' % (BDM1+bubble)-P0 mixed FEM
            [u,p,w,edge,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option);
    end
    % compute error
    tic;
    % ================== error for velocity ==================
    if isfield(pde,'exactu') 
        switch elemType
            case 'RT0-P0'
                erruL2(k) = getL2errorRT0(node,elem,pde.exactu,u);
                % interpolation
                uI = faceinterpolate(pde.exactu,node,edge,'RT0');
                ufreeDof = eqn.ufreeDof;
                u0 = u(ufreeDof);
                uI0 = uI(ufreeDof);
                erruIuhH1(k) = sqrt((u0-uI0)'*eqn.A*(u0-uI0));
                erruInf(k) = max(abs(u0-uI0));
            case 'BDM1B-P0'
                erruL2(k) = getL2errorBDM1(node,elem,pde.exactu,u);
                % interpolation
                uI = u;
                uI(1:2*size(edge,1)) = faceinterpolate(pde.exactu,node,edge,'BDM1');
                ufreeDof = eqn.ufreeDof;
                u0 = u(ufreeDof);
                uI0 = uI(ufreeDof);
                erruIuhH1(k) = sqrt((u0-uI0)'*eqn.A*(u0-uI0));
                erruInf(k) = max(abs(u0-uI0));
        end
    end
    % ================== error for pressure ==================
    if isfield(pde,'exactp')
        area = simplexvolume(node,elem);
        errpL2(k) = getL2error(node,elem,pde.exactp,p);
        pI = Lagrangeinterpolate(pde.exactp,node,elem,'P0');
        rp = recoverP02P1(node,elem,p,'LA');
        % --------------------------
        errpL2(k) = getL2error(node,elem,pde.exactp,p);
        errpL2re(k) = getL2error(node,elem,pde.exactp,rp);
        errpIphL2(k) = sqrt(dot((pI-p).^2,area));
        errpInf(k) = max(abs(p-pI));
    end
    % ================== error for vorticity ==================
    if isfield(pde,'exactw')
        switch elemType
            case 'RT0-P0'
                errwL2(k) = getL2error(node,elem,pde.exactw,w);
                wI = Lagrangeinterpolate(pde.exactw,node,elem,'P1');
                errwIwh(k) = sqrt((w-wI)'*eqn.Mv*(w-wI));
            case 'BDM1B-P0'
                errwL2(k) = getL2error(node,elem,pde.exactw,w);
                wI = Lagrangeinterpolate(pde.exactw,node,elem,'P2',edge);
                errwIwh(k) = sqrt((w-wI)'*eqn.Mv*(w-wI));
        end
    end
    % ================== record ==================
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = size(node,1);
    if option.plotflag && N(k) < 3e3 % show mesh and solution for small size
       figure(1); showresult(node,elem,p); pause(0.1);
       figure(2); 
       showsolutionRT(node,elem,u); pause(0.1);
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    tic;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    meshTime(k) = toc;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.80,0.40]);
    subplot(1,3,1)
    h = 1./sqrt(N(1:k));
    showrateh3(h,erruL2(1:k),2,'k-+','||u-u_h||',...
               h,erruIuhH1(1:k),2,'r-*','||u_I-u_h||_1',...
               h,erruInf(1:k),2,'b-*','||u_I-u_h||_{\infty}');
    subplot(1,3,2)
    showrateh4(h,errpL2(1:k),2,'k-+', '||p - p_h||',...
               h,errpL2re(1:k),2,'g-+','||p- p_h^r||',...
               h,errpIphL2(1:k),2,'r-+','||p_I - p_h||',...
               h,errpInf(1:k),2,'b-+','||p_I - p_h||_{\infty}');
    subplot(1,3,3)
    showrateh2(h,errwL2(1:k),2,'k-+','||w - w_h||',...
               h,errwIwh(1:k),2,'r-+','||w_I - w_h||');
end

% Output
err = struct('N',N,'uL2',erruL2(1:k),'uInf',erruInf,'uIuhH1',erruIuhH1(1:k),...
             'pL2',errpL2(1:k),'pIphL2',errpIphL2(1:k),'pInf',errpInf(1:k),'pL2re',errpL2re(1:k),...
             'wL2',errwL2(1:k),'wIwhL2',errwIwh(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assmble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error on screen
ts = zeros(k,3); ts = char(ts);
% error of u
fprintf('===========================================================\n');
display('   N   ||u_I-u_h||_1  ||u-u_h||   ||u_I-u_h||_{max}');
display([num2str(err.N) ts num2str(err.uIuhH1,'%0.5e') ts num2str(err.uL2,'%0.5e') ts num2str(err.uInf,'%0.5e')]);
% error of p
fprintf('===========================================================\n');
display('   N   ||p_I-p_h||   ||p-p_h||    ||p_I-p_h||_{max}      ||p_I - p^r_h||');
display([num2str(err.N) ts num2str(err.pIphL2,'%0.5e') ts num2str(err.pL2,'%0.5e') ts num2str(err.pInf,'%0.5e') ts num2str(err.pL2re,'%0.5e')]);
% error of w
fprintf('===========================================================\n');
display('   N   ||w_I-w_h||   ||w-w_h||');
display([num2str(err.N) ts num2str(err.wIwhL2,'%0.5e') ts num2str(err.wL2,'%0.5e')]);
fprintf('===========================================================\n');
