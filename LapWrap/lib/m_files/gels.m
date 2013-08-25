function [x] = gels(a,b,param)
%
% GELS   Solves an overdetermined or underdetermined linear problem.
%
%    [X] = GELS(A,B) same as [X] = GELS(A,B,'default') 
%
%    [X] = GELS(A,B,PARAM) if A is m-by-n and m > n then GELS computes the solution of the
%    linear least squares problem
%            min || A * X - B ||_F ,
%    if A is m-by-n and m < n then GELS computes the solution of the minimum norm problem
%            min || X ||_F s.t. A*X = B
%    PARAM controls the algorithm used:
%	
%       'dc'
%             performs SVD factorization of A using the divide and conquer
%             method then solve the undetermined or overdetermined problem.
%             calls LAPACK routine _GELSD
%
%       'svd'
%             performs SVD factorization of A using QR iterations then solve the
%             undetermined or overdetermined problem.
%             calls LAPACK routine _GELSS
%
%       'y'
%             performs QR factorization of A then solve the undetermined or
%             overdetermined problem.
%             calls LAPACK routine _GELSY
%
%       'default'
%             (default)
%             performs QR factorization of A without pivoting, then solve the
%             undetermined or overdetermined problem.  The fact that there is no
%             pivoting in the QR factorization means that A(:,1:min(size(A))
%             needs to be full rank.
%             calls LAPACK routine _GELS
%
%    When A is m-by-n with m > n then [X] = GELS(A,B) is the same as X=A\B.
%
%    When A is m-by-n with m < n  then [X] = GELS(A,B,PARAM) is _not_ the same as
%    X=A\B. [X] = GELS(A,B,PARAM) solves the problem
%    (*)    min_X || X ||_F s.t. AX=B
%    while X=A\B solves AX=B by first performing a QR factorization of A
%          [Q,R,E] = QR(A)
%    then X is obtained via
%          X = E(:,1:m)*(R(:,1:m)\(Q'*b))] 
%	The \ method is faster than GELS, X is such that AX=B but X is not (in the
%	general case) the solution of (*).
%
%    See also: \, MLDIVIDE.
%
	x=[];
	

	if (nargout~=1),
		disp('Wrong number of output parameters');
		return;
	end;

	if (nargin~=2 && nargin~=3),
		disp('Wrong number of input parameters');
		return;
	end;
	
	%Do all the boring checking
	if (~isnumeric(a)),
		disp('The matrix is not composed of numeric values.');
		return;
	elseif (isinteger(a)),
		a=double(a);
	end;

	if ((nargin == 2)||(strcmp(param,'default'))),
		trans = 'N';
		[x,info] = gels_2(trans,a,b);
	else
		if (strcmp(param,'dc')),
			rcond = -1; % machine precision
			[x,rank,info] = gelsd(a,b,rcond);
		elseif (strcmp(param,'svd')),
			rcond = -1; % machine precision
			[x,rank,info] = gelss(a,b,rcond);
		elseif (strcmp(param,'y')),
			rcond = eps; % machine precision
			[x,rank,info] = gelsy(a,b,rcond);
		else
			disp('Type of algorithm not correct. Please choose dc, svd, y, or default.');
			return;
		end;
	end;	

	%Do post-call boring checking
	if (info~=0),
		disp('Resolution failed, no solution returned!');
		x=[];
		return;
	end;
