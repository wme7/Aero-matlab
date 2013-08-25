function [m,val,vec,isuppz,info] = syevr(jobz,range,uplo,a,rl,ru,abstol)
%SYEV   Computing eigen information.
%   [M,VAL,VEC,ISUPPZ,INFO] = SYEVR(JOBZ,RANGE,UPLO,A,RL,RU,ABSTOL) 
%   computes eigenvalues and, optionally, eigenvectors of a symmetric matrix A.
%
%   JOBZ allows the user to choose if eigenvectors are to be computed.
%   JOBZ
%          = 'N':  Compute eigenvalues only;
%          = 'V':  Compute eigenvalues and eigenvectors.
%
%   RANGE 
%          = 'A': all eigenvalues will be found.
%          = 'V': all eigenvalues in the half-open interval (RL,RU]
%                 will be found.
%          = 'I': the RL-th through RU-th eigenvalues will be found.
%
%   UPLO allows the user to choose which part of the matrix will be referenced.
%   UPLO
%          = 'U':  Upper triangle of A is stored;
%          = 'L':  Lower triangle of A is stored.
%
%   RL and RU T
%       The upper and lower bound for searching eigenvalues.
%       They are integer or real, depending on the value of the RANGE parameter.
%
%   M 
%       The number of computed eigenvalues.
%
%   ISUPPZ
%
%       The support of the eigenvectors in Z, i.e., the indices
%       indicating the nonzero elements in Z. The i-th eigenvector
%       is nonzero only in elements ISUPPZ( 2*i-1 ) through
%       ISUPPZ( 2*i ).
%       Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
%
%   ABSTOL 
%       The absolute error tolerance for the eigenvalues.
%
%   If only eigenvalues are desired, the QR algorithm is used. 
%
%   No input parameter are optional. If you want an easy to use interface,
%   see SYEV_DRIVER.

    vec = [];
	val = [];
	info = [];
	m =	[];
	isuppz = [];

	if (nargin~=7),
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
	ldz=int32(n);
	
	if (range=='A'),
		vl=0;vu=0;il=int32(0);iu=int32(0);
	elseif (range=='V'),
		vl=rl;vu=ru;il=int32(0);iu=int32(0);
	elseif (range=='I'),
		il=rl;iu=ru;vl=int32(0);vu=int32(0);
	else
		disp('Bad value for parameter RANGE.');
		return;
	end;
	
	if (isComplex(a)),
		if (isDouble(a))
			abstol=-1;
			w=zeros(double(n),1);
			work=complex(zeros(1));
			rwork=zeros(1,1);
			iwork=int32(zeros(1,1));
			lrwork=int32(-1);
			liwork=int32(-1);
			lwork=int32(-1);
			info=int32(-1);
			isuppz=int32(zeros(double(2*max(1,n)),1));
			z=complex(zeros(double(ldz),double(max(1,n))));
			m=int32(zeros(1,1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info]=lapack_zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
			
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			lrwork=int32(rwork(1));

			work=complex(zeros(double(lwork),1));
			iwork=int32(zeros(double(lwork),1));
			rwork=zeros(double(lwork),1);
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info]=lapack_zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
			
		else
			abstol=single(-1);
			w=single(zeros(double(n),1));
			work=single(complex(zeros(1)));
			rwork=single(zeros(1,1));
			iwork=int32(zeros(1,1));
			lrwork=int32(-1);
			liwork=int32(-1);
			lwork=int32(-1);
			info=int32(-1);
			isuppz=int32(zeros(double(2*max(1,n)),1));
			z=single(complex(zeros(double(ldz),double(max(1,n)))));
			m=int32(zeros(1,1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info]=lapack_cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
			
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			lrwork=int32(rwork(1));

			work=single(complex(zeros(double(lwork),1)));
			iwork=int32(zeros(double(liwork),1));
			rwork=single(zeros(double(lrwork),1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info]=lapack_cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);

		end;
	else
		if (isDouble(a)),
			abstol=-1;
			w=zeros(double(n),1);
			work=zeros(1);
			iwork=int32(zeros(1));
			lwork=int32(-1);
			liwork=int32(-1);
			info=int32(-1);
			isuppz=int32(zeros(double(2*double(max(1,n))),1));
			z=zeros(double(ldz),double(max(1,n)));
			m=int32(zeros(1,1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork,  iwork, liwork, info]=lapack_dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			work=zeros(double(lwork),1);
			iwork=int32(zeros(double(liwork),1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork,  iwork, liwork, info]=lapack_dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
		else
			abstol=single(-1);
			w=single(zeros(double(n),1));
			work=single(zeros(1));
			iwork=int32(zeros(1));
			lwork=int32(-1);
			liwork=int32(-1);
			info=int32(-1);
			isuppz=int32(zeros(double(2*max(1,n)),1));
			z=single(zeros(double(ldz),double(max(1,n))));
			m=int32(zeros(1,1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork,  iwork, liwork, info]=lapack_ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
			lwork=int32(work(1));
			liwork=int32(iwork(1));
			work=single(zeros(double(lwork),1));
			iwork=int32(zeros(double(liwork),1));
			[jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, val, vec, ldz, isuppz, work, lwork,  iwork, liwork, info]=lapack_ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
		end;
	end;
