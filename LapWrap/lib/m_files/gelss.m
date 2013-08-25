function [x,rank,info] = gelss(a,b_in,rcond)
%GELSD   Computing eigen information.
%   [X,INFO] = GELSS(A,B,RCOND) 
%
%	No input parameter are optional. If you want an easy to use interface,
%	see GELS.

    x = [];
	rank = [] ;
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
	nrhs=int32(size(b_in,2));
	lda=int32(m);
	ldb=int32(max(m,n));

	if (isDouble(b_in)),
		b=zeros(double(ldb),double(nrhs));
		b(1:size(b_in,1) , 1:size(b_in,2)) = b_in;
	else
		b=single(zeros(double(ldb),double(nrhs)));
		b(1:size(b_in,1) , 1:size(b_in,2)) = b_in;
	end;
	rank=int32(0);

	if (isComplex(a) || isComplex(b)),
		if (isDouble(a) || isDouble(b))

			if (~isComplex(b)),
				b=complex(b);
			end;
			if (~isComplex(a)),
				a=complex(a);
			end;
			if (~isDouble(b)),
				b=double(b);
			end;
			if (~isDouble(a)),
				a=double(a);
			end;
			s=zeros(double(min(m,n)),1);
			work=complex(zeros(1));
			lrwork=int32(5*min(m,n));
			rwork=zeros(double(lrwork),1);
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, info]=lapack_zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
			lwork=int32(work(1));
			work=complex(zeros(double(lwork),1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, info]=lapack_zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);

		else %complex single
			if (~isComplex(b)),
				b=complex(b);
			end;
			if (~isComplex(a)),
				a=complex(a);
			end;
			s=single(zeros(double(min(m,n)),1));

			work=single(complex(zeros(1,1)));
			work(1)=sqrt(-1)+1;
			lwork=int32(single(-1));
			lrwork=int32(5*min(m,n));
			rwork=single(zeros(double(lrwork),1));
			info=int32(single(-1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, info]=lapack_cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
			lwork=int32(work(1));
			work=single(complex(zeros(double(lwork),1)));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, info]=lapack_cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);

		end;
	else %real double
		if (isDouble(a) || isDouble(b)),
			if (~isDouble(b)),
				b=double(b);
			end;
			if (~isDouble(a)),
				a=double(a);
			end;
			s=zeros(double(min(m,n)),1);
			
			work=zeros(1);
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, info]=lapack_dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
			lwork=int32(work(1));
			work=zeros(double(lwork),1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, info]=lapack_dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
		else %real single
			s=single(zeros(double(min(m,n)),1));

			work=single(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, info]=lapack_sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
			lwork=int32(work(1));
			work=single(zeros(double(lwork),1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, info]=lapack_sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
		end;
	end;

	x=x(1:double(n),1:double(nrhs));

