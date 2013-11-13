function lp = lagrange_poly ( xi, i )

%LAGRANGE_POLY   determine the Lagrange polynomial associated with the 
%                i-th interpolating point of a given set of interpolating
%                points
%
%     calling sequences:
%             lp = lagrange_poly ( xi, i )
%             lagrange_poly ( xi, i )
%
%     inputs:
%             xi      vector containing the interpolating points
%             i       index of the interpolating point with which the
%                     Lagrange polynomial is to be associated
%
%     output:
%             lp      vector containing coefficients of the Lagrange
%                     polynomial associated with the i-th interpolating 
%                     point of the given set of interpolating points
%
%     NOTE:
%             this routine is primarily a support routine used by LAGRANGE
%             to determine the Lagrange form of the interpolating polynomial
%
%             to evaluate the Lagrange polynomial, apply the MATLAB function
%             'polyval' to the output from this routine 
%

n = length ( xi );
temp = [1];
denom = 1;

if ( ( i < 1 ) | ( i > n ) )
   disp ( 'lagrange_poly error: index i out of bounds' )
   return
elseif ( i == 1 )
   for j = 2:n
       temp = conv ( temp, [1 -xi(j)] );
	   denom = denom * ( xi(1) - xi(j) );
   end
elseif ( i == n )
   for j = 1:n-1
       temp = conv ( temp, [1 -xi(j)] );
	   denom = denom * ( xi(n) - xi(j) );
   end
else
   for j = 1:i-1
       temp = conv ( temp, [1 -xi(j)] );
	   denom = denom * ( xi(i) - xi(j) );
   end
   for j = i+1:n
       temp = conv ( temp, [1 -xi(j)] );
	   denom = denom * ( xi(i) - xi(j) );
   end
end

if ( nargout == 0 )
   disp(temp/denom)
else
   lp = temp/denom;
end