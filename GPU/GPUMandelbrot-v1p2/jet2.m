function cmap = jet2(m)
% Jet colormap with added fade to black

%   Copyright 2010 The Mathworks, Inc.
%   $Revision: 1$
%   $Date: 2010-11-08$

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

% A list of break-point colors
colors = [
    0.0  0.0  0.5
    0.0  0.0  1.0
    0.0  0.5  1.0
    0.0  1.0  1.0
    0.5  1.0  0.5
    1.0  1.0  0.0
    1.0  0.5  0.0
    1.0  0.0  0.0
    0.5  0.0  0.0    
    0.5  0.0  0.0
    1.0  0.0  0.0
    1.0  0.5  0.0
    1.0  1.0  0.0
    0.5  1.0  0.5
    0.0  1.0  1.0
    0.0  0.5  1.0
    0.0  0.0  1.0
    0.0  0.0  0.5
    0.0  0.0  0.0
    ];

% Now work out the indices into the map
N = size( colors, 1 );
idxIn = 1:N;
idxOut = linspace( 1, N, m );
cmap = [
    interp1( idxIn, colors(:,1), idxOut )
    interp1( idxIn, colors(:,2), idxOut )
    interp1( idxIn, colors(:,3), idxOut )
    ]';