function [val,vec,info] = syev(jobz,uplo,a)
%SYEV   Computing eigen information.
%   [VAL,VEC,INFO] = SYEV(JOBZ,UPLO,A) computes eigenvalues and,
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
			lrwork=int32(max(1,3*n-2));
			rwork=zeros(double(lrwork),1);
			lwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,info]=lapack_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork,info);
			lwork=int32(work(1));
			work=complex(zeros(double(lwork),1));
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,info]=lapack_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork,info);
		else
			w=single(zeros(double(n),1));
			work=single(complex(zeros(1,1)));
			work(1)=sqrt(-1)+1;
			lrwork=int32(single(max(1,3*n-2)));
			rwork=single(zeros(double(lrwork),1));
			lwork=int32(single(-1));
			info=int32(single(-1));
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,info]=lapack_cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork,info);
			lwork=int32(work(1));
			work=single(complex(zeros(double(lwork),1)));
			[jobz, uplo, n, vec, lda, val, work, lwork, rwork,info]=lapack_cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork,info);

		end;
	else
		if (isDouble(a)),
			w=zeros(double(n),1);
			work=zeros(1);
			lwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, info]=lapack_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
			lwork=int32(work(1));
			work=zeros(double(lwork),1);
			[jobz, uplo, n, vec, lda, val, work, lwork, info]=lapack_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
		else
			w=single(zeros(double(n),1));
			work=single(zeros(1));
			lwork=int32(-1);
			info=int32(-1);
			[jobz, uplo, n, vec, lda, val, work, lwork, info]=lapack_ssyev(jobz, uplo, n, a, lda, w, work, lwork, info);
			lwork=int32(work(1));
			work=single(zeros(double(lwork),1));
			[jobz, uplo, n, vec, lda, val, work, lwork, info]=lapack_ssyev(jobz, uplo, n, a, lda, w, work, lwork, info);
		end;
	end;
