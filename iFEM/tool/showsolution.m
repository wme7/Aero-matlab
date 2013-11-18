function showsolution(node,elem,u,varargin)
%% SHOWSOLUTION plots the solution u on a triangular mesh in 2-D.
%
%    showsolution(node,elem,u) displays the functoin u on a topological
%    2-dimensional mesh given by node and elem matrices. The function u
%    could be piecewise constant or piecewise linear. 
%
%    showsolution(node,elem,u,viewangle) changes the display angle. The
%    deault view angle on planar meshes is view(2) and view(3) for surface
%    meshes. 
%
%    showsolution(node,elem,u,'param','value','param','value'...) allows
%    additional patch param/value pairs to be used when displaying the
%    mesh. 
%
%   Example:
%     f = inline('sin(2*pi*x).*cos(2*pi*y)');
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%         [node,elem] = uniformrefine(node,elem);
%     end
%     u = f(node(:,1),node(:,2));
%     subplot(1,3,1);
%     showsolution(node,elem,u);
%     subplot(1,3,2);
%     showsolution(node,elem,u,[-62,58]);
%     subplot(1,3,3);
%     showsolution(node,elem,u,'EdgeColor','k');
%
%   See also showmesh, showsolution3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if size(u,1)< size(u,2), u = u'; end

if (length(u) == size(elem,1)) % piecewise constant functoin
%     x = node(:,1);
%     y = node(:,2);
%     h = patch(x(elem'), y(elem'), u(elem'), u(elem'),'FaceColor','interp','EdgeColor', 'interp');
    NT = size(elem,1);
    elemnew = reshape(1:3*NT,NT,3);
    nodenew = node(elem(:),:);
    unew = repmat(u,3,1);
    h = trisurf(elemnew, nodenew(:,1), nodenew(:,2), unew', ...
            'FaceColor', 'interp', 'EdgeColor', 'interp');
    axis equal; axis tight;
end

if (length(u) == size(node,1)) % piecewise linear functoin
    dim = size(node,2);
    nv = size(elem,2);
    if (dim==2) && (nv==3)    % planar triangulation
        h = trisurf(elem, node(:,1), node(:,2), u', ...
                'FaceColor', 'interp', 'EdgeColor', 'interp');
%         view(2);
    end
    if (dim==2) && (nv==4)    % planar quadrilateration
        x = node(:,1);
        y = node(:,2);
        h = patch(x(elem'), y(elem'), u(elem'), u(elem'),'FaceColor', ...
                  'interp', 'EdgeColor', 'interp');
%         axis equal; axis tight;
    end    
    if (dim==3)     % surface mesh
        h = trisurf(elem,node(:,1),node(:,2),node(:,3),u','FaceColor',...
                    'interp', 'EdgeColor', 'interp');
        view(3);
    end
end

if nargin>3
    if isnumeric(varargin{1})
        view(varargin{1});
        if nargin>4
            set(h,varargin{2:end});
        end
    else
        set(h,varargin{1:end});        
    end
end
% axis equal; axis tight;