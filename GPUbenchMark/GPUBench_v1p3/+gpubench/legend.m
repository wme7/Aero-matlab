function l = legend( varargin )
%legend  Create a legend and make it match the background

%   Copyright 2011 The MathWorks, Inc.

l = legend( varargin{:} );
figh = ancestor( l, 'figure' );
set( l, 'Color', get( figh, 'Color' ) );
