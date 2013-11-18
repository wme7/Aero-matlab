function option = mgoptions(option,N)
%% MGOPTIONS default options for multigrid solver
%
% x0 = 0;	tol = 1e-8;   solver = 'CG';    preconditioner = 'Vcycle';
% N0 = 500; coarsegridsolver = 'direct';    mu = 1; 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~isfield(option,'x0')
    option.x0 = zeros(N,1);    
end 
if ~isfield(option,'tol')
    option.tol = 1e-8;
end 
if ~isfield(option,'solvermaxit')
    option.solvermaxit = min(N,200);
end
if isfield(option,'solver')
    option.solver = upper(option.solver);
else
    option.solver = 'CG';
end
if isfield(option,'smoother')
    option.smoother = upper(option.smoother);
else
    option.smoother = 'GS';  % defacult one is Gauss-Seidel
    % option.solver = 'JAC' % Jacobi preconditioner
end
if ~isfield(option,'smoothingstep')  % smoothing steps
    option.smoothingstep = 1;
end
if ~isfield(option,'smoothingparameter')  % smoothing parameter
    option.smoothingparameter = 1;
end
if isfield(option,'preconditioner')
    option.preconditioner = upper(option.preconditioner);
elseif ~strcmp(option.solver,'V') && ~strcmp(option.solver,'W') && ~strcmp(option.solver,'F')
    option.preconditioner = 'V';
end
if ~isfield(option,'N0')
    option.N0 = 50;
end
if ~isfield(option,'coarsegridsolver')
    option.coarsegridsolver = 'direct';    % solver in the coarsest grid
end
if ~isfield(option,'printlevel')
    option.printlevel = 1;
end
if ~isfield(option,'setupflag') 
    option.setupflag = true;
end