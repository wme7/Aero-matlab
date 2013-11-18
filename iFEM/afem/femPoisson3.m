function [err,time,solver,eqn] = femPoisson3(node,elem,pde,bdFlag,option,varargin)
%% FEMPOISSON3 solve Poisson equation by various finite element methods
%
%   FEMPOISSON3 computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if nargin >=1 && ischar(node)
    option.elemType = node;
    clear node
end
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = cubemesh([0,1,0,1,0,1],0.25); % default mesh is a cube
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = sincosdata3;                          % default data
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
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'P1'     % piecewise linear function P1 element
            [u,Du,eqn,info] = Poisson3(node,elem,pde,bdFlag,option);
        case 'CR'     % piecewise linear function CR element
            [u,Du,eqn,info] = Poisson3CR(node,elem,pde,bdFlag,option);
        case 'P2'     % piecewise quadratic function
            [u,Du,eqn,info] = Poisson3P2(node,elem,pde,bdFlag,option);
        case 'WG'     % weak Galerkin element
            [u,Du,eqn,info] = Poisson3WG(node,elem,pde,bdFlag,option);            
    end
    % compute error
    tic;
    if isfield(pde,'Du')
        if ~isempty(Du)
            errH1(k) = getH1error3(node,elem,pde.Du,Du);
        else
            errH1(k) = getH1error3(node,elem,pde.Du,u);            
        end
    end
    if isfield(pde,'exactu')
        errL2(k) = getL2error3(node,elem,pde.exactu,u);        
        % interpolation
        switch elemType
            case 'P1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'CR'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.face);
            case 'P2'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
            case 'WG'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.face);
        end
        erruIuh(k) = sqrt((u-uI)'*eqn.A*(u-uI));
        errMax(k) = max(abs(u-uI));
    end
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(u);
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
    showrate2(N(1:k),errH1(1:k),1,'-*','||Du-Du_h||',...
              N(1:k),errL2(1:k),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:k),erruIuh(1:k),1,'m-+','||Du_I-Du_h||',...
              N(1:k),errMax(1:k),1,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N,'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assmble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));

%% Display error
ts = zeros(k,3); ts = char(ts);
display(' #Dof   ||u-u_h||     ||Du-Du_h||   ||DuI-Du_h||  ||uI-u_h||_{max}');
display([num2str(err.N) ts num2str(err.L2,'%0.5e') ts num2str(err.H1,'%0.5e')...
         ts num2str(err.uIuhH1,'%0.5e') ts num2str(err.uIuhMax,'%0.5e')]);