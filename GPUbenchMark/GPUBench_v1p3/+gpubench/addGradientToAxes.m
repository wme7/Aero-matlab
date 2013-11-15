function addGradientToAxes( axh )
%addGradientToAxes  Add some shading to the axes
%
%   gpubench.addGradientToAxes(axh)

%   Copyright 2011 The MathWorks, Inc.

% Pick colors for the top and bottom of the canvas
topCol = [0.87 0.87 0.87];
botCol = [0.98 0.98 0.98];

% First, remove any existing shading
p = findall( axh, 'Tag', 'FunkyAxes:AxesBackground' );
if ~isempty( p )
    delete( p );
end
xlim = get( axh, 'XLim' );
ylim = get( axh, 'YLim' );
if strcmpi( get( axh, 'YDir' ), 'Reverse' )
    ylim = fliplr( ylim );
end
zlim = get( axh, 'ZLim' );
clim = get( axh, 'CLim' );
% Make sure the grid appears over the patch
set( axh, 'Layer', 'top' );

% Now create the patch
r = [ topCol(1)*[1 1]  botCol(1)*[1 1] ];
g = [ topCol(2)*[1 1]  botCol(2)*[1 1] ];
b = [ topCol(3)*[1 1]  botCol(3)*[1 1] ];
cdata = cat( 3, r, g, b );
% A bug in patch means that it sometimes resets CLim!
ptch = patch( ...
    'Parent', axh, ...
    'XData', [xlim,fliplr(xlim)], ...
    'YData', [ylim(1)*[1 1] ylim(2)*[1 1]], ...
    'ZData', -inf*ones(1,4), ...
    'CData', cdata, ...
    'FaceColor', 'interp', ...
    'FaceLighting', 'none', ...
    'Tag', 'FunkyAxes:AxesBackground', ...
    'HandleVisibility', 'off', ...
    'HitTest', 'off' );
set( axh, 'CLim', clim );
% Make sure the patch is at the bottom!
iSendToBack( axh, ptch )
hold( axh, 'on' );
% Listen for changes in the limits so that we don't see patch edges
hgax = handle( axh );
setappdata( axh, 'LimitListeners', {
    handle.listener( hgax, findprop( hgax, 'XLim' ), 'PropertyPostSet', @iUpdatePatch )
    handle.listener( hgax, findprop( hgax, 'YLim' ), 'PropertyPostSet', @iUpdatePatch )
    } );


    function iUpdatePatch(~,~)
        if ~isempty( ptch ) && ishandle( ptch )
            xlim = get( axh, 'XLim' );
            ylim = get( axh, 'YLim' );
            zlim = get( axh, 'ZLim' );
            set( ptch, ...
                'XData', [xlim,fliplr(xlim)], ...
                'YData', [ylim(1)*[1 1] ylim(2)*[1 1]], ...
                'ZData', zlim(1)*ones(1,4) );
        end
    end % iUpdatePatch

end % addGradientToAxes

function iSendToBack( axh, bg )
% Send an HG object to the back of the display stack
children = findall( axh,'Parent',axh );
set( axh, 'Children', [children(children~=bg);bg] )
end % iSendToBack