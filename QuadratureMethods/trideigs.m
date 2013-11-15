function [v,lambda] = trideigs ( a, b, TOL, Nmax, vects )

%QRST           determine all of the eigenvalues (and optionally all of the
%               eigenvectors) of a symmetric tridiagonal matrix using the 
%               QR algorithm with Wilkinson shift
%
%     calling sequences:
%             [lambda, v] = qrst ( a, b, TOL, Nmax, vects )
%             [lambda, v] = qrst ( a, b, TOL, Nmax )
%             lambda = qrst ( a, b, TOL, Nmax )
%             qrst ( a, b, TOL, Nmax, vects )
%             qrst ( a, b, TOL, Nmax )
%
%     inputs:
%             a       vector containing elements along the main diagonal
%                     of the symmetric tridiagonal matrix whose eigenvalues
%                     are to be determined
%             b       vector containing elements along the off diagonal
%                     of the symmetric tridiagonal matrix whose eigenvalues
%                     are to be determined
%             TOL     convergence tolerance
%             Nmax    maximum number of iterations
%             vects   optional input argument
%                     matrix containing eigenvector information produced
%                     during the reduction of the original symmetric
%                     matrix to symmetric tridiagonal form
%                     - this input is needed only if computation of the 
%                       eigenvectors is requested (by including the second
%                       output argument) and the original matrix was not in
%                       symmetric tridiagonal form
%
%     output:
%             lambda  vector containing the eigenvalues of the symmetric
%                     tridiagonal matrix determined by the vectors a and b
%             v       optional output argument
%                     matrix containing the eigenvectors of the symmetric
%                     tridiagonal matrix determined by the vectors a and b
%                     - the i-th column of this matrix is an eigenvector
%                       which corresponds to the i-th eigenvalue in the
%                       vector lambda
%                     - eigenvectors will be not computed if this second
%                       output argument is omitted
%
%     NOTE:
%             if the maximum number of iterations is exceeded, a message
%             to this effect will be displayed, along with the number of
%             eigenvalues which had been determined - these eigenvalues
%             will be returned in the last entries of the output vector
%             lambda
%

n = length(a);
if ( length(b) == n-1 ) b(2:n) = b(1:n-1); end;

c = zeros ( 1, n );
s = zeros ( 1, n );
shift = 0;
togo = n;

if ( nargout >= 2 )
   if ( nargin >= 5 )
      v = vects;
   else
      v = eye(n);
   end;
end;

for its = 1 : Nmax

    if ( togo == 1 )
	   lambda(1) = a(1) + shift;
	   disp ( its );
	   return;
	end;
	
    trace = a(togo-1) + a(togo);
	det   = a(togo-1)*a(togo) - b(togo)*b(togo);
	disc  = sqrt ( trace*trace - 4*det );
	mu1 = (1/2) * ( trace + disc );
	mu2 = (1/2) * ( trace - disc );
	if ( abs ( mu1 - a(togo) ) < abs ( mu2 - a(togo) ) )
	   s = mu1;
	else
	   s = mu2;
	end;

    shift = shift + s;
	for i = 1:togo 
	    a(i) = a(i) - s;
	end;
	
	oldb = b(2);
    for i = 2:togo
        j = i-1;
	    r = sqrt ( a(j)^2 + oldb^2 );
	    c(i) = a(j) / r;
	    s(i) = oldb / r;
	    a(j) = r;
	    temp1 = c(i)*b(i) + s(i)*a(i);
	    temp2 = -s(i)*b(i) + c(i)*a(i);
	    b(i) = temp1;
	    a(i) = temp2;
	    if ( i ~= togo ) oldb = b(i+1); b(i+1) = c(i)*b(i+1); end;
    end;

    a(1) = c(2)*a(1) + s(2)*b(2);
    b(2) = s(2)*a(2);
    for i = 2:togo-1
        a(i) = s(i+1)*b(i+1) + c(i)*c(i+1)*a(i);
	    b(i+1) = s(i+1)*a(i+1);
    end;
    a(togo) = c(togo)*a(togo);
	
	if ( nargout >= 2 )
	   for i = 2 : togo
	       col1 = v(:,i-1) * c(i) + v(:,i) * s(i);
		   v(:,i) = -s(i) * v(:,i-1) + c(i) * v(:,i);
		   v(:,i-1) = col1;
	   end;
	end;
	
	if ( abs(b(togo)) < TOL )
	   lambda(togo) = a(togo) + shift;
	   disp([lambda(togo) its]);
	   togo = togo - 1;
	end;

end;

disp ( 'qrst error: Maximum number of iterations exceeded' );
disp ( sprintf ( '%d eigenvalues determined \n', n-togo ) );