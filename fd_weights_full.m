%-----------------------------------------------------------------------
      function[c] = fd_weights_full(xx,x,m);
%
%     This routine evaluates the derivative based on all points
%     in the stencils.  It is more memory efficient than "fd_weights"
%
%     This set of routines comes from the appendix of 
%     A Practical Guide to Pseudospectral Methods, B. Fornberg
%     Cambridge Univ. Press, 1996.   (pff)
%
%     Input parameters:
%       xx -- point at wich the approximations are to be accurate
%       x  -- array of x-ordinates:   x(0:n)
%       m  -- highest order of derivative to be approximated at xx
%
%     Output:
%       c  -- set of coefficients c(0:n,0:m).
%             c(j,k) is to be applied at x(j) when
%             the kth derivative is approximated by a 
%             stencil extending over x(0),x(1),...x(n).
%
%
%     UPDATED 8/26/03 to account for matlab "+1" index shift.
%     Follows p. 168--169 of Fornberg's book.
%
%
      n1 = length(x);
      n  = n1-1;
      m1 = m+1;

      c1       = 1.;
      c4       = x(1) - xx;

      c = zeros(n1,m1);
      c(1,1) = 1.;

      for i=1:n;                   i1  = i+1;
         mn = min(i,m);            mn1 = mn+1;
         c2 = 1.;
         c5 = c4;
         c4 = x(i1)-xx;
         for j=0:i-1;              j1  = j+1;
            c3 = x(i1)-x(j1);
            c2 = c2*c3;
            for k=mn:-1:1;         k1  = k+1;
               c(i1,k1) = c1*(k*c(i1-1,k1-1)-c5*c(i1-1,k1))/c2;
            end;
            c(i1,1) = -c1*c5*c(i1-1,1)/c2;
            for k=mn:-1:1;         k1  = k+1;
               c(j1,k1) = (c4*c(j1,k1)-k*c(j1,k1-1))/c3;
            end;
            c(j1,1) = c4*c(j1,1)/c3;
         end;
         c1 = c2;
      end;