function showsolutionCR(node,elem,u,varargin)

NT = size(elem,1);

%% Generate a big triangulation
elemnew = reshape(1:3*NT,NT,3);
nodenew = node(elem(:),:);

%% Evaluate piecewise linear function
elem2edge = dofedge(elem);
unew(1:NT) = -u(elem2edge(:,1)) + u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(NT+(1:NT)) = u(elem2edge(:,1)) - u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(2*NT+(1:NT)) = u(elem2edge(:,1)) + u(elem2edge(:,2)) - u(elem2edge(:,3));

%% Plot the solution
showsolution(nodenew,elemnew,unew,varargin{:});