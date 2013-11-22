function m = minmod(a)
%% minmod function
% In the following implementation, it is assumed 'a' is a vector with 3
% elements, 1<=j<=l where l = 3.
% by Manuel Diaz, NTU, 2013.11.08

signtest = sign(a(1,:)) == sign(a(2,:)) & sign(a(1,:)) == sign(a(3,:));
    s = sign(a(1,:)); 
    m = s.*min(abs(a)).*signtest;