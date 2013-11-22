function m_tilde = MODminmod(a,M,h)
%% MODminmod function
% In the following implementation, it is assumed 'a' is a vector with 3
% elements, 1<=j<=l where l = 3.
% by Manuel Diaz, NTU, 2013.11.08

if abs(a(1,:)) <= M*h^2
    m_tilde = a(1,:);
else
    m_tilde = minmod(a);
end