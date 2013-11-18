function [elemType,maxIt,maxN,L0,refType,option] = mfemoption(option)

if isfield(option,'elemType')
    elemType = upper(option.elemType);
else
    elemType = 'RT0';    
end

if isfield(option,'maxIt')
    maxIt = option.maxIt;
else
    maxIt = 4;    
end

if isfield(option,'maxN')
    maxN = option.maxN;
else
    maxN = 1e5;    
end

if isfield(option,'L0')
    L0 = option.L0;
else
    L0 = 0;    
end

if isfield(option,'refType')
    refType = lower(option.refType);
else
    refType = 'red';    
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