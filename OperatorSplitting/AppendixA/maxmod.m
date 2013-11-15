function s=maxmod(a,b)
% The maxmod limiter
s=0.5*(sign(a)+sign(b)).*max(abs(a),abs(b));