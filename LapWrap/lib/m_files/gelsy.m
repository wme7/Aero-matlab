function [x,rank,info] = gelsy(a,b_in,rcond)
%GELSY   Computing eigen information.
%   [X,RANK,INFO] = GELSY(A,B,RCOND) 
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

	jpvt=int32(zeros(double(n),1));

	rank=int32(-1);

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
			work=complex(zeros(1));
			lrwork=int32(2*n);
			rwork=zeros(double(lrwork),1);
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, rwork, info]=lapack_zgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
			lwork=int32(work(1));
			work=complex(zeros(double(lwork),1));
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, rwork, info]=lapack_zgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);

		else %complex single
			if (~isComplex(b)),
				b=complex(b);
			end;
			if (~isComplex(a)),
				a=complex(a);
			end;

			work=single(complex(zeros(1,1)));
			work(1)=sqrt(-1)+1;
			lwork=int32(single(-1));
			lrwork=int32(2*n);
			rwork=single(zeros(double(lrwork),1));
			info=int32(single(-1));
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, rwork, info]=lapack_cgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
			lwork=int32(work(1));
			work=single(complex(zeros(double(lwork),1)));
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, rwork, info]=lapack_cgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);

		end;
	else %real double
		if (isDouble(a) || isDouble(b)),
			if (~isDouble(b)),
				b=double(b);
			end;
			if (~isDouble(a)),
				a=double(a);
			end;

			work=zeros(1);
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a_dest, lda, x, ldb, jpvt, rcond, rank, work, lwork, info]=lapack_dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
			lwork=int32(work(1));
			work=zeros(double(lwork),1);
			[ m, n, nrhs, a_dest, lda, x, ldb, jpvt, rcond, rank, work, lwork, info]=lapack_dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
		else %real single

			work=single(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, info]=lapack_sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
			lwork=int32(work(1));
			work=single(zeros(double(lwork),1));
			[ m, n, nrhs, a, lda, x, ldb, jpvt, rcond, rank, work, lwork, info]=lapack_sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
		end;
	end;

    %Do post-call boring checking
    if (info~=0),
        disp('Resolution failed, no solution returned!');
        x=[];
        return;
    end;


	x=x(1:double(n),1:double(nrhs));

