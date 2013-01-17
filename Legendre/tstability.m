function tstability(p, form)
  
  % p= order of Jameson-Schmidt-Turkell Runge-Kutta scheme
  % form = format to plot stability curve(s) in Matlab
  

  % create polynomial coefficients for contraction factor
  mypoly = zeros(p+1,1);
  for j=0:p-1
    mypoly(j+1) = (1/factorial(p-j));
  end 
  
  
  % plot dependence of polynomial roots on theta 
  for ntheta = 1:100
    theta = 2*pi*(ntheta-1)/(100-1);
    
    mypoly(p+1) = 1+ exp(i*theta);
    
    ro      = roots(mypoly);
    
    hold on; plot(real(ro), imag(ro), form); hold off;
    xlabel('real(alpha)'); ylabel('imag(alpha)');
  end
