function [x,rank,info] = gelsd(a,b_in,rcond)
%GELSD   Computing eigen information.
%   [X,INFO] = GELSD(A,B,RCOND) 
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
			rwork=zeros(1);
			iwork=int32(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, iwork, info]=lapack_zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
			lwork=int32(work(1));
			work=complex(zeros(double(lwork),1));
			lrwork=int32(rwork(1));
			smlsiz=25;
			nlvl = max( 0, int32( log2( double(min( m,n ))/double(smlsiz+1) ) ) + 1 );
			lrwork = 10*min(m,n) + 2*min(m,n)*smlsiz + 8*min(m,n)*nlvl + 3*smlsiz*nrhs + (smlsiz+1)*(smlsiz+1);
		
			liwork = max(1, 3*min(m,n)*nlvl + 11*min(m,n));
			rwork=zeros(double(lrwork),1);
			iwork=int32(zeros(double(liwork),1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, iwork, info]=lapack_zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);

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
			rwork=single(zeros(1));
			iwork=int32(zeros(1));
			info=int32(single(-1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, iwork, info]=lapack_cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
			lwork=int32(work(1));
			work=single(complex(zeros(double(lwork),1)));
			lrwork=int32(rwork(1));
			smlsiz=25;
			nlvl = max( 0, int32( log2( double(min( m,n ))/double(smlsiz+1) ) ) + 1 );
			lrwork = 10*min(m,n) + 2*min(m,n)*smlsiz + 8*min(m,n)*nlvl + 3*smlsiz*nrhs + (smlsiz+1)*(smlsiz+1);
		
			liwork = max(1, 3*min(m,n)*nlvl + 11*min(m,n));
			rwork=single(zeros(double(lrwork),1));
			iwork=int32(zeros(double(liwork),1));
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, rwork, iwork, info]=lapack_cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);

		end;
	else %real double
		if (isDouble(a) || isDouble(b)),
			smlsiz=int32(32); 
			nlvl = max( 0, int32( log2( double(min( m,n ))/double((double(smlsiz)+1)) )  + 1) );
			lwork= int32(11*double(max(m,n)) + 2*double(max(m,n))*double(smlsiz) + 8*double(max(m,n))*double(nlvl) + double(max(m,n))*double(nrhs));
			liwork=int32(3 * double(min(m,n)) * double(nlvl) + 11 * double(min(m,n)));
			iwork=int32(zeros(double(liwork),1));
			if (~isDouble(b)),
				b=double(b);
			end;
			if (~isDouble(a)),
				a=double(a);
			end;
			s=zeros(double(min(m,n)),1);
			info=int32(-1);
			iwork=int32(zeros(double(liwork),1));
			work=zeros(double(lwork),1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, iwork, info]=lapack_dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
		else %real single
			smlsiz=int32(25);
			nlvl = max( 0, int32( log2( double(min( m,n ))/double((smlsiz+1)) ) ) + 1 );
			liwork=int32(3 * min(m,n) * nlvl + 11 * min(m,n));
			lwork= int32(11*max(m,n) + 2*max(m,n)*smlsiz + 8*max(m,n)*nlvl + max(m,n)*nrhs);
			iwork=int32(zeros(double(liwork),1));
			work=single(zeros(double(lwork),1));
			s=single(zeros(double(min(m,n)),1));

			info=int32(-1);
			[ m, n, nrhs, a, lda, x, ldb, s, rcond, rank, work, lwork, iwork, info]=lapack_sgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
		end;
	end;

	x=x(1:double(n),1:double(nrhs));

