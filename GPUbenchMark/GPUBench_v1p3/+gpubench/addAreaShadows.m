function addAreaShadows( axh )
%addAreaShadows  add shadows to any patches in the axes
%
%   gpubench.addAreaShadows(axh)

%   Copyright 2011 The MathWorks, Inc.

% The color to use for the shadow. Should probably adapt to the background,
% but for now let's just hard-code it.
shadowColor = 0.6*[1 1 1];

% The shadow size is expressed as a ratio of the axes width/height
shadowSize = 1/300;

% First remove any existing ones
l = findall( axh, 'Tag', 'FunkyAxes:AreaShadow' );
if ~isempty( l )
    delete( l )
end

% Work out the axes limits and correspopnding shadow offsets
dx = diff( get( axh, 'XLim' ) )*shadowSize;
dy = -diff( get( axh, 'YLim' ) )*shadowSize;
if strcmpi( get( axh, 'XDir' ), 'Reverse' )
    dx = -dx;
end
if strcmpi( get( axh, 'YDir' ), 'Reverse' )
    dy = -dy;
end
z0 = get( axh, 'ZLim' );
z0 = z0(1) - max(0.1,diff( z0 )*shadowSize);

% Find all existing patches in the axes, excluding any background
patches = findobj( axh, 'type', 'patch' );
patches( strcmpi( get( patches, 'Tag' ), 'FunkyAxes:AxesBackground' ) ) = [];
props_to_copy = {
    'XData'
    'YData'
    'ZData'
    'FaceColor'
    'EdgeColor'
    'LineStyle'
    'LineWidth'
    };

% Now go through each patch, copying various properties into a new "shadow"
% patch.
vals = cell( size( props_to_copy ) );
l = -1*ones( size( patches ) );
for jj=1:numel( patches )
    l(jj) = patch( nan, nan, nan, ...
        'Tag', 'FunkyAxes:AreaShadow', ...
        'Parent', axh );
    for kk=1:numel( props_to_copy )
        vals{kk} =  get( patches(jj), props_to_copy{kk} );
    end
    args = [props_to_copy,vals]';
    set( l(jj), args{:} );
    
    set( l(jj), ...
        'FaceColor', shadowColor, ...
        'EdgeColor', shadowColor, ...
        'XData', get( patches(jj), 'XData' ) + dx, ...
        'YData', get( patches(jj), 'YData' ) + dy, ...
        'ZData', z0*ones( size( get( patches(jj), 'XData' ) ) ), ...
        'HitTest', 'off' );
end
