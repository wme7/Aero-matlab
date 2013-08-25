function [val,vec,info] = syevd(jobz,uplo,a)
%SYEV   Computing eigen information.
%   [VAL,VEC,INFO] = SYEVD(JOBZ,UPLO,A) computes eigenvalues and,
%	optionally, eigenvectors of a symmetric matrix A.
%
%	JOBZ allows the user to choose if eigenvectors are to be computed.
%	JOBZ
%          = 'N':  Compute eigenvalues only;
%          = 'V':  Compute eigenvalues and eigenvectors.
%
%	UPLO allows the user to choose which part of the matrix will be referenced.
%	UPLO
%          = 'U':  Upper triangle of A is stored;
%          = 'L':  Lower triangle of A is stored.
%
%   If only eigenvalues are desired, the QR algorithm is used. 
%
%	No input parameter are optional. If you want an easy to use interface,
%	see SYEV_DRIVER.

    vec = [];
	val = [];
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

    if (size(a,1)~=size(a,2)),
        disp('The matrix must be square.');
        return;
	end;
	n=int32(size(a,1));
	lda=int32(n);
	
	if (isComplex(a)),
		if (isDouble(a))

			w=zeros(double(n),1);
			work=complex(zeros(1));
			rwork=zeros(1,1);
			iwork=int32(zeros(1,1));
			lrwork=int32(-1);
			liwork=int32(-1);
			lwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,lrwork, iwork, liwork, info]=lapack_zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork,lrwork, iwork, liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			lrwork=int32(rwork(1));

			work=complex(zeros(double(lwork),1));
			iwork=int32(zeros(double(lwork),1));
			rwork=zeros(double(lwork),1);
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,lrwork, iwork, liwork, info]=lapack_zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork,lrwork, iwork, liwork, info);
		else
			w=single(zeros(double(n),1));
			work=single(complex(zeros(1)));
			rwork=single(zeros(1,1));
			iwork=int32(zeros(1,1));
			lrwork=int32(-1);
			liwork=int32(-1);
			lwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,lrwork, iwork, liwork, info]=lapack_cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork,lrwork, iwork, liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			lrwork=int32(rwork(1));

			work=single(complex(zeros(double(lwork),1))); 
			iwork=int32(zeros(double(liwork),1));
			rwork=single(zeros(double(lrwork),1));
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,lrwork, iwork, liwork, info]=lapack_cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork,lrwork, iwork, liwork, info);

		end;
	else
		if (isDouble(a)),
			w=zeros(double(n),1);
			work=zeros(1);
			iwork=int32(zeros(1));
			lwork=int32(-1);
			liwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, iwork,liwork, info]= lapack_dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork,liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			work=zeros(double(lwork),1);
			iwork=int32(zeros(double(liwork),1));
			[jobz, uplo, n, vec, lda, val, work, lwork, iwork,liwork, info]= lapack_dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork,liwork, info);
		else
			w=single(zeros(double(n),1));
			work=single(zeros(1));
			iwork=int32(zeros(1));
			lwork=int32(-1);
			liwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, iwork,liwork, info]= lapack_ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork,liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			work=single(zeros(double(lwork),1));
			iwork=int32(zeros(double(liwork),1));
			[jobz, uplo, n, vec, lda, val, work, lwork, iwork,liwork, info]= lapack_ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork,liwork, info);
		end;
	end;
