 function[M1,M2,M3] = integrate_vf(k,w,f,v)
      M1   = (k*sum(v .* w .* f));    % Number density
      M2   = (k*sum((v .^2 ) .* w .* f));   % Macrospic moment in x
      M3   = (k*sum(1/2*( v .*abs(v).^2 ).* w .* f)); % Energy Density
 end