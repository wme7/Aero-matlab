function option = femoption(option)
% [elemType,maxIt,maxN,L0,refType,option] = femoption(option)

if ~isfield(option,'elemType')
    option.elemType = 'P1';    
end

if ~isfield(option,'maxIt')
    option.maxIt = 4;    
end

if ~isfield(option,'maxN')
    option.maxN = 1e5;    
end

if ~isfield(option,'L0')
    option.L0 = 0;    
end

if ~isfield(option,'refType')
    option.refType = 'red';
end

if ~isfield(option,'printlevel')
    option.printlevel = 1;    
end

if ~isfield(option,'plotflag')
    option.plotflag = 1;    
end

if ~isfield(option,'rateflag')
    option.rateflag = 1;    
end

if ~isfield(option,'tol')
    option.tol = 1e-8;    
end

if ~isfield(option,'lumpflag') % mass lumping
    option.lumpflag = 0;       % no mass lumping
end