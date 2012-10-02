function s=minmod(a,b)
% The minmod limiter
s=0.5*(sign(a)+sign(b)).*min(abs(a),abs(b));