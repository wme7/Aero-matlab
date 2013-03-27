function [XIter,XHistory,FXHistory,failureFlag] = Newton(X,F,X0,varargin)
%% Newton's method for solving a system of nonlinear equations
%
% Newton(X,F,X0) solves nonlinear system F(X)=0 by Newton's method, using
% the given initial approximation X0.
%
% [XIter,XHistory,FXHistory,failureFlag] = Newton(X,F,X0,opts);
% Required inputs are X, F, and X0. Optional inputs are:
%   opts.maxIter,  the max number of iterations (default 5);
%   opts.tol, tolerance (default 1e-15),
%             the method stops if max(abs(F(X))) < opts.tol.
% The required output is XIter, which is the approximate solution
% Optional outputs are:
%   XHistory, the history of all computed approximate solutions;
%   FXHistory simply is F(XHistory);
%   failureFlag is 0 if max(abs(F(XIter))) < opts.tol, for all finite
%               components of F(XIter), or 1 otherwise.
%
% Class support for inputs X and F sym, X0 double, opts structure.
% Class support for all outputs: double
%
%%   Examples:
%
% disp('Scalar example:');
% clear all; X=sym('X'); F=sin(X); X0=1;
% Sol=Newton(X,F,X0) % is equivalent to
% Sol=Newton(X,F,X0,struct('tol',1e-15,'maxIter',5)) % also equivalent to
% opts.tol=1e-15; opts.maxIter=5; Sol=Newton(X,F,X0,opts)
%
% disp('The same example, only using different notaion');
% clear all; t=sym('t'); f=sin(t); t0=1;
% Sol=Newton(t,f,t0)
%
% disp('More equations (2) than unknowns (1), system consistent');
% clear all; X=sym('X'); F=[sin(X); tan(X)]; X0=1;
% [Sol,XHistory,FXHistory]=Newton(X,F,X0)
%
% disp('More equations (2) than unknowns (1), system inconsistent');
% clear all; X=sym('X'); F=[sin(X); tan(X)-1]; X0=1;
% [Sol,XHistory,FXHistory]=Newton(X,F,X0)
%
% disp('2 equations and 2 unknowns, system consistent');
% clear all; X=sym('X',[2 1]); F=[sin(X(1)); tan(X(2))-1]; X0=ones(2,1);
% [Sol,XHistory,FXHistory]=Newton(X,F,X0)
%
% disp('2 equations and 2 unknowns, system consistent');
% clear all; X=sym('X',[2 1]); X0=2*ones(2,1);
% F=[6*X(1)^3+X(1)*X(2)-3*X(2)^3-4; ...
%    X(1)^2-18*X(1)*X(2)^2+16*X(2)^3+1];
% [Sol,XHistory,FXHistory]=Newton(X,F,X0)

% License: BSD
% Copyright (c) 2011-2012 A.V. Knyazev, Andrew.Knyazev@ucdenver.edu
% http://math.ucdenver.edu/~aknyazev/
% $Revision: 1.2 $  $Date: 6-Jan-2012
% Tested in MATLAB 7.13 (R2011b) and its Symbolic Math Toolbox 5.7 (R2011b)
% Does NOT work in Octave 3.4.2 and below because of poor symbolic toolbox

failureFlag = 1;
%check the required input parameters for consistency
if nargin < 3,
    error('Newton:NotEnoughInputs',...
        strcat('There must be at least 3 input agruments: ',...
        ' X, F and X0.'));
end%if
if ~isa(X,'sym'),
    fprintf('The first required input parameter, X, is class %s.\n',...
        class(X));
    error('Newton:NotSymbolicInputX',...
        'The first required input parameter, X, must be class sym.');
end%if
if ~isa(F,'sym'),
    fprintf('The second required input parameter, F, is class %s.\n',...
        class(F));
    error('Newton:NotSymbolicInputF',...
        'The second required input parameter, F, must be class sym.');
end%if
if ~all(isfinite(X0)),
    fprintf('The third required input parameter, X0, is class %s.\n',...
        class(X0));
    display(X0);
    error('Newton:NotFiniteInputX0',...
        'The third required input parameter, X0, must be finite double.');
end%if
if size(X,2)~=1,
    fprintf('The first required input parameter, X,  has the size');
    disp(size(X));
    error('Newton:NotColumnVectorX',...
        strcat('The first required input parameter, X, must be a scalar',...
        ' or a column vector.'));
end%if
if size(X0,2)~=1,
    fprintf('The third required input parameter, X0,  has the size');
    disp(size(X0));
    error('Newton:NotColumnVectorX0',...
        strcat('The third required input parameter, X0, must be a scalar',...
        ' or a column vector.'));
end%if
X0 = cast(full(X0),'double'); % cast X0 into full double
if ~all(size(X)==size(X0)),
    fprintf('The first required input parameter, X,  has the size.');
    disp(size(X));
    fprintf('The third required input parameter, X0, has the size.');
    disp(size(X0));
    error('Newton:SizeMismatchXandX0',...
        strcat('The first and the third required input parameters, ',...
        'X and X0, must be of the same size.'));
end%if

maxIter=5; tol=1e-15; % setting the defaults
if nargin == 4, % overwting the defaults, if input is present
    opts=varargin{1};
    if isfield(opts,'maxIter'),
        if ~isfloat(opts.maxIter) || ~isscalar(opts.maxIter) || ...
                ~isreal(opts.maxIter) || (opts.maxIter<0.5),
            fprintf('The optional parameter opts.maxIter is ');
            disp(opts.maxIter);
            error('Newton:InvalidOptsmaxIter',...
                strcat('The optional parameter opts.maxIter, ',...
                'must be a positive scalar integer.'));
        end%if
        maxIter=round(opts.maxIter);
        if maxIter ~= opts.maxIter
            fprintf('The parameter opts.maxIter %g is rounded to %g .\n',...
                opts.maxIter, maxIter);
        end%if
    end%if
    if isfield(opts,'tol')
        if ~isfloat(opts.tol) || ~isscalar(opts.tol) || ...
                ~isreal(opts.tol) || (opts.tol<=0),
            fprintf('The optional parameter opts.tol is ');
            disp(opts.tol);
            error('Newton:InvalidOptsTol',...
                strcat('The optional parameter opts.tol, ',...
                'must be a positive scalar float.'));
        end%if
        tol = cast(full(opts.tol),'double');
    end%if
end%if

H = jacobian(F,X);
XIter=X0; FXIter=subs(F,X,XIter);
XHistory=XIter; FXHistory=FXIter;
if  ~all(isfinite(FXIter)),
    warning('Newton:NotFiniteFunctionValueAtX0',...
        'All function values at X0 must be finite. Iterations stopped.');
    return
end%if

if max(abs(FXIter)) < tol,
    failureFlag = 0; % required tolerance already, no iterations
    XHistory=XIter; FXHistory=FXIter;
    return
end%if
    
% preallocation for Hisory
XHistory=zeros(size(X,1),maxIter+1);
FXHistory=zeros(size(F,1),maxIter+1);

for j=1:maxIter, % the main iterative loop
    XHistory(:,j)=XIter; FXHistory(:,j)=FXIter;
    HIter=double(subs(H,X,XIter));
    if  ~all(isfinite(HIter)),
        warning('Newton:NotFiniteJacobianValue',...
            'All Jacobian values must be finite. Iterations stopped.');
        break;
    end%if
    XIter=XIter-HIter\FXIter;
    if  ~all(isfinite(XIter)),
        warning('Newton:NotFiniteXIterValue',...
            'All values of XIter must be finite. Iterations stopped.');
        break;
    end%if
    FXIter=double(subs(F,X,XIter));
    if  ~all(isfinite(FXIter))
        warning('Newton:NotFiniteFunctionValue',...
            strcat('All function values FXIter at XIter ',....
            'must be finite. Iterations stopped.'));
        break;
    end%if
    if max(abs(FXIter)) < tol,
        failureFlag = 0; break;
    end%if
    XIterDiff=XIter-XHistory(:,j);
    if  max(abs((XIterDiff))) < eps,
        warning('Newton:StagnationDetected',...
            'Stagnation, no progress in the approximate solution');
        break
    end%if
end%for

XHistory(:,j+1)=XIter; FXHistory(:,j+1)=FXIter; % to store the last step
XHistory(:,j+2:maxIter+1)=[]; FXHistory(:,j+2:maxIter+1)=[]; %remove preallocated
end%function