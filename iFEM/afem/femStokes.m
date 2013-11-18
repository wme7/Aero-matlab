function [err,time,solver,eqn] = femStokes(node,elem,pde,bdFlag,option,varargin)
%% FEMSTOKES solve the Stokes equation by various finite element methods
%
%   FEMSTOKES computes approximations to the Stokes equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% Created by Ming Wang at Nov., 2012.
%
% See also femPoisson femrateStokes
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Default setting of mesh and pde data
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = squaremesh([0,1,0,1],0.125);  % default mesh is a square
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = Stokesdata1;                          % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end

%% Parameters
option = femoption(option);
maxIt = option.maxIt;   maxN = option.maxN; L0 = option.L0;
option = femStokesoption(option);
elemType = option.elemType; refType = option.refType;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(option.refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(option.refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
erruH1 = zeros(maxIt,1); errpL2 = zeros(maxIt,1); 
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch upper(elemType)
        case 'P2P1'
            [u,p,edge,A,eqn,info] = StokesP2P1(node,elem,pde,bdFlag,option);
        case 'P2P0'
            [u,p,edge,A,eqn,info] = StokesP2P0(node,elem,pde,bdFlag,option);
        case 'ISOP2P1'
            [u,p,edge,A,eqn,info] = StokesisoP2P1(node,elem,pde,bdFlag,option);
        case 'ISOP2P0'
            [u,p,edge,A,eqn,info] = StokesisoP2P0(node,elem,pde,bdFlag,option);
        case 'CRP0'
            [u,p,edge,A,eqn,info] = StokesCRP0(node,elem,pde,bdFlag,option);
        case 'CRP1'
            [u,p,edge,A,eqn,info] = StokesCRP1(node,elem,pde,bdFlag,option);
        case 'MINI'
            [u,p,edge,A,eqn,info] = StokesMini(node,elem,pde,bdFlag,option);
        case 'P1BP1'
            [u,p,edge,A,eqn,info] = StokesP1bP1(node,elem,pde,bdFlag,option);
    end    
    % compute error
    tic;
%     if strcmp(elemType,'P1BP1')
%         uI = Lagrangeinterpolate(pde.exactu,node,elem);
%     elseif strcmp(elemType(1:2),'CR')
%         uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',edge);        
%     end
    if strcmp(elemType,'P2P0') || strcmp(elemType,'P2P1') || ...
       strcmp(elemType,'isoP2P0') || strcmp(elemType,'isoP2P1')
        uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    elseif strcmp(elemType,'CRP0') || strcmp(elemType,'CRP1')
        uI = pde.exactu((node(edge(:,1),:)+node(edge(:,2),:))/2);
    elseif strcmp(elemType,'Mini')
        uI = pde.exactu(node);
    elseif strcmp(elemType,'P1bP1')
        % bubble part won't be taken in the error computation
        Nv = size(node,1); NT = size(elem,1);
        u0 = pde.exactu(node);
        uI = u;
        uI([(1:Nv)'; NT+Nv+(1:Nv)']) = u0(:);
    end
    erruH1(k) = sqrt((u-uI(:))'*A*(u-uI(:)));
    errpL2(k) = getL2error(node,elem,pde.exactp,p);
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),erruH1(k), errpL2(k));
    end
    % record solver information
    itStep(k) = info.itStep;
%    stopErr(k) = info.stopErr;
%    flag(k) = info.flag;
    % plot
    N(k) = size(node,1);
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
        if length(p) == size(elem,1) % piecewise constant function
            p = recoverP02P1(node,elem,p);
        end
        figure(1);  showresult(node,elem,p);
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
h = 1./sqrt(N(1:k));

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh(h,erruH1(1:k),1,'k-+','|u_I-u_h|_1');
    subplot(1,2,2)
    showrateh(h,errpL2(1:k),1,'m-+','||p-p_h||');
end

%% Output
err = struct('N',N,'uH1',erruH1(1:k),'pL2',errpL2(1:k));          
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assmble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error
ts = zeros(k,6); ts = char(ts);
display('#nodes     |u_I-u_h|_1      ||p-p_h||  ');
% format shorte
display([num2str(err.N) ts num2str(err.uH1,'%0.5e') ts num2str(err.pL2,'%0.5e')]);