function lp = lagrange ( xi, fi )

%LAGRANGE     determine the Lagrange form of the interpolating polynomial 
%             associated with a given set of interpolating points and
%             function values
%
%     calling sequences:
%             lp = lagrange ( xi, fi )
%             lagrange ( xi, fi )
%
%     inputs:
%             xi      vector containing the interpolating points
%             fi      vector containing function values
%                     the i-th entry in this vector is the function
%                     value associated with the i-th entry in the 'xi'
%                     vector
%
%     output:
%             lp      vector containing coefficients of the Lagrange
%                     form of the interpolating polynomial associated with 
%                     the given set of interpolating points and function
%                     values
%
%     NOTE:
%             to evaluate the Lagrange form, apply the MATLAB function
%             'polyval' to the output from this routine 
%

n = length ( xi );
m = length ( fi );

if ( n ~= m )
   disp ( 'lagrange error: number of ordinates and number of function values must be equal' )
   return
end

temp = [0];
for i = 1:n
    temp = temp + lagrange_poly ( xi, i ) * fi(i);
end

if ( nargout == 0 )
   disp(temp)
else
   lp = temp;
end