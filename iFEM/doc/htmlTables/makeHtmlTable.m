function sOut = makeHtmlTable(M, T, rowNames, colHeaders, colors, strPrecision)
%MAKEHTMLTABLE  Display matrix contents as an HTML table
%   makeHtmlTable(M) where M is an array
%   makeHtmlTable(M, T) where T is a cell array of strings equal in size
%   to M. The contents of T are displayed preferentially when T is not empty.
%   makeHtmlTable(M, T, rowNames, colHeaders) where colHeaders and
%   rowNames are cell arrays.
%   makeHtmlTable(M, T, rowNames, colHeaders, colors) where colors is a
%   standard three-column colormap. Color is scaled like in IMAGESC. NaNs
%   in M are mapped to the color white.
%   makeHtmlTable(M, T, rowNames, colHeaders, colors, strPrecision) where
%   strPrecision specifies the digits of precision displayed as explained
%   in MAT2STR.
%
%   Example:
%   makeHtmlTable([1 2; 3 4])

%   Author: Ned Gulley
%   Modified by Long Chen.
%   Copyright 2009 The MathWorks, Inc.

if nargin < 6
    % Use full precision
    strPrecision = 15;
end

[nr,nc] = size(M);
TM = cell(size(M));
for i = 1:nr
    for j = 1:nc
        TM{i,j} = mat2str(M(i,j),strPrecision);
    end
end

if (nargin < 2) || isempty(T)
    T = TM;
else
    [nrt,nct] = size(T);
    if (nrt ~= nr) || (nct ~= nc)
        error('Input matrices M and T must be the same size')
    end
    for i = 1:nr
        for j = 1:nc
            if isempty(T{i,j})
                T{i,j} = TM{i,j};
            end
        end
    end
end

if nargin < 4
    colHeaders = [];
    rowNames = [];
end

if (nargin < 5) || isempty(colors)
    colorFlag = false;
else
    colorFlag = true;
end


% I threw this color flag in for fun.
if colorFlag
    minM = nanmin(M(:));
    maxM = nanmax(M(:));
    if maxM == minM
        colorFlag = false;
    else
        normM = (M-minM)/(maxM-minM);
        [nColors,three] = size(colors);
        if three ~= 3
            error('Colormap must have three columns')
        end
        % Turn normM into an integer map into the color table
        normM = floor(normM*(nColors-1))+1;

        % Any cells marked with NaN will be white;
        normM(isnan(normM)) = nColors + 1;
        colors = [colors; 1 1 1];
    end
end


s = {};

s{end+1} = sprintf('<html><table border="1" cellpadding="4" cellspacing="0">\n');

if ~isempty(colHeaders)
    s{end+1} = sprintf('<tr>');
    if ~isempty(rowNames)
        s{end+1} = sprintf('<td></td>');
    end
    for j = 1:nc
        s{end+1} = sprintf('<td>');
        s{end+1} = sprintf('%s',colHeaders{j});
        s{end+1} = sprintf('</td>');
    end
    s{end+1} = sprintf('</tr>\n');
end


for i = 1:nr
    s{end+1} = sprintf('<tr>');
    if ~isempty(rowNames) % add by Long Chen
        s{end+1} = sprintf('<td>%s</td>',rowNames{i});
    end

    for j = 1:nc
        if colorFlag
            s{end+1} = sprintf('<td bgcolor=#%s%s%s>', ...
                dec2hex(floor(255*colors(normM(i,j),1)),2), ...
                dec2hex(floor(255*colors(normM(i,j),2)),2), ...
                dec2hex(floor(255*colors(normM(i,j),3)),2));
        else
            s{end+1} = sprintf('<td>');
        end
        s{end+1} = sprintf('%s', T{i,j});
        s{end+1} = sprintf('</td>');
    end
    s{end+1} = sprintf('</tr>\n');
end

s{end+1} = sprintf('</table></html>\n');

if nargout==1
    sOut = s;
else
    for i = 1:length(s)
        fprintf('%s',s{i});
    end
end


% =========================================================================
function hexStr = rgb2hexStr(color)
hexStr = dec2hex(floor(255*color),2);

