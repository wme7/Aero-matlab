function [x,info] = gels_2(trans,a,b)
%GELS_2   Computing eigen information.
%   [X,INFO] = GELS_2(TRANS,A,B) solves overdetermined or underdetermined
%	real linear systems involving an M-by-N matrix A, or its transpose,
%	using a QR or LQ factorization of A.  It is assumed that A has full
%	rank.
%
%	The following options are provided:
%
%	1. If TRANS = 'N' and m >= n:  find the least squares solution of
%	   an overdetermined system, i.e., solve the least squares problem
%	   minimize || B - A*X ||.
%
%	2. If TRANS = 'N' and m < n:  find the minimum norm solution of
%	   an underdetermined system A * X = B.
%
%	3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
%	   an undetermined system A**T * X = B.
%
%	No input parameter are optional. If you want an easy to use interface,
%	see GELS.

    x = [];
	info = [];
	
	if (nargin~=3),
		disp('Wrong number of input parameters.');
		return;
	end;
	
	    %Do all the boring checking
	if (~isnumeric(a)),
		disp('The matrix is not composed of numeric values.');
		return;
	elseif (isinteger(a)),
		a=double(a);
	end;

	m=int32(size(a,1));
	n=int32(size(a,2));
	nrhs=int32(size(b,2));
	lda=int32(m);
	ldb=int32(max(m,n));
	
	if (isComplex(a)),
		if (isDouble(a))

			b_local=complex(zeros(double(ldb),double(nrhs)));
			b_local(1:size(b,1), 1:size(b,2))=double(b);
			work=complex(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_zgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
			lwork=int32(work(1));
			work=complex(zeros(double(lwork),1));
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_zgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);

		else %complex single
			b_local=complex(single(zeros(double(ldb),double(nrhs))));
			b_local(1:size(b,1), 1:size(b,2))=single(b);
			work=single(complex(zeros(1,1)));
			work(1)=sqrt(-1)+1;
			lwork=int32(single(-1));
			info=int32(single(-1));
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_cgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
			lwork=int32(work(1));
			work=single(complex(zeros(double(lwork),1)));
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_cgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);

		end;
	else %real double
		if (isDouble(a)),
			b_local=zeros(double(ldb),double(nrhs));
			b_local(1:size(b,1), 1:size(b,2))=double(b);
			work=zeros(1);
			lwork=int32(-1);
			info=int32(-1);
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_dgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
			lwork=int32(work(1));
			work=zeros(double(lwork),1);
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_dgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
		else %real single
			b_local=single(zeros(double(ldb),double(nrhs)));
			b_local(1:size(b,1), 1:size(b,2))=single(b);
			work=single(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_sgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
			lwork=int32(work(1));
			work=single(zeros(double(lwork),1));
			[ trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info]=lapack_sgels(trans, m, n, nrhs, a, lda, b_local, ldb, work, lwork,info);
		end;
	end;

	clear x;
	
	if (strcmp(trans,'N'))
		x=b(1:double(n),1:double(nrhs));
	else
		x=b(1:double(m),1:double(nrhs));
	end;
