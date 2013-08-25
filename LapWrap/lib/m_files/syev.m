function [vec,val] = syev(a,param)
%
% SYEV   Eigenvalues and eigenvectors of symmetric/Hermitian matrix.
%
%	[D] = SYEV(A) computes the eigenvalues of a symmetric/Hermitian matrix A.
%	First A is reduced to tridiagonal form so that A=Q1*T*Q1' then the
%	symmetric QR iteration algorithm is used on the tridiagonal matrix T to
%	find D, the eigenvalues of T. D corresponds to the eigenvalues of A. This
%	routine calls LAPACK routine _SYEV/_HEEV.
%
%	[V,D] = SYEV(A) computes eigenvalues and eigenvectors of a
%	symmetric/Hermitian matrix A. First A is reduced to tridiagonal form so
%	that A = Q1*T*Q1' then the symmetric QR iteration algorithm is used on the
%	tridiagonal matrix T to find D, the eigenvalues of T, and Q2 the
%	eigenvectors of T so that T = Q2*D*Q2'. D corresponds to the eigenvalues of
%	A. The eigenvectors of A are such that V = Q1*Q2.
%	
%	[V,D] = SYEV(A,PARAM) computes eigenvalues and eigenvectors of a
%	symmetric/Hermitian matrix A. Same as [V,D] = SYEV(A) except that PARAM
%	enables to choose the method used in the symmetric tridiagonal
%	eigenproblem.
%	PARAM controls the algorithm used:
%
%       'dc'
%             stands for "Divide and Conquer"
%             calls LAPACK routine _SYEVD/_HEEVD
%
%       'mrrr'
%             stands for "Multiple Relatively Robust Representations"
%             calls LAPACK routine _SYEVR/_HEEVR
%
%       'bx'
%             stands for "bisection and inverse iteration"
%             calls LAPACK routine _SYEVX/_HEEVX
%
%       'qr'
%             (default)
%             stands for " symmetric QR iteration"
%             calls LAPACK routine _SYEV/_HEEV
%
%    When A is symmetric/Hermitian SYEV(A,'qr') is the same as EIG(A).
%
%    See also: EIG.
%
	vec=[];
	val=[];
	
	uplo = 'U';

	if (nargout<=1),
		jobz='N';
	elseif (nargout==2),
		jobz='V';
	else
		disp('Wrong number of output parameters');
		return;
	end;

	%Do all the boring checking
	if (~isnumeric(a)),
		disp('The matrix is not composed of numeric values.');
		return;
	elseif (isinteger(a)),
		a=double(a);
	end;
	
	if (size(a,1)~=size(a,2)),
		disp('The matrix must be square.');
		return;
	end;

	if ((nargin == 1)||(nargout<=1)||(strcmp(param,'qr'))),
		[val,vec,info] = syev_2(jobz,uplo,a);
	else
		if (strcmp(param,'dc')),
			[val,vec,info] = syevd(jobz,uplo,a);
		elseif (strcmp(param,'bx')),
			abstol=-1;
			range='A';
			rl = 0; ru = 0;
			[m,val,vec,ifail,info] = syevx(jobz,range,uplo,a,rl,ru,abstol);
		elseif (strcmp(param,'mrrr')),
			abstol=-1;
			range='A';
			rl = 0; ru = 0;
			[m,val,vec,isuppz,info] = syevr(jobz,range,uplo,a,rl,ru,abstol);
		else
			disp('Type of algorithm not correct. Please choose dc, bii or rrr.');
			return;
		end;
	end;	

	%Do post-call boring checking
	if (info~=0),
		disp('Resolution failed, no eigen information returned!');
		vec=[];
		val=[];
		return;
	end;

	if (nargout==2)
		val=diag(val);
	end
	if (nargout==1)
		val=vec;
	end
