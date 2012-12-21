% Plot u0
%if plot_figs == 1; plot(x,u0); axis([0,2,0,1.2]); end;
if plot_figs == 1; plot(x,u0); axis('tight'); end;

% Define coeficients matrix
a = diag([1 12 180 2800 44100 698544 11099088 176679360]);

% transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
ut = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    for j = 1:nx
    i = l+1;            % Dummy index
    P = LegMat(k,xi);   % Legendre Matrix
    ut(i,j) = a(i,i).*sum(w(:,j).*u0(:,j).*P(:,i));
    end
end

%% transform degress of freedom u(t)_{l,i} back to space-time values u(x,t)
u = (ut'*P')';
figure
plot(x,u); axis('tight');