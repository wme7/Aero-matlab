%% Printing Variables to HTML Tables in Published Code
% One often-requested feature for Cell Mode publishing is to have a command
% that simply displays the output of a MATLAB® variable as a table. This
% file shows how to do that.
%
% Copyright 2009 The MathWorks, Inc.

%% Passing through HTML
% Beginning with R2007a, you can pass HTML directly through to your output
% file by wrapping it in an HTML tag.
%
% Look at the difference between these two lines to see how it works.

disp('Option 1: This is <b>bold</b>')

%%

disp('<html>Option 2: This is <b>bold</b></html>')

%%
% Here's a matrix...

m = magic(4)

%%
% which I can now display with the help of a function like makeHtmlTable.

makeHtmlTable(m);

%% 
% Here is the help entry for makeHtmlTable

help makeHtmlTable

%% Add column headers and row names
% If you have them, we can show them.

r = {'row 1','row 2','row 3','row 4'};
c = {'col 1','col 2','col 3','col 4'};
makeHtmlTable(m,[],r,c);

%% Colorize the table cells
% Here's a little hack for matrix visualization.

colorFlag = true;
cmap = gray;
cmap(:,1) = 1;
makeHtmlTable(m,[],r,c,cmap);

%% Some test cases
% If all these run to completion, we're doing pretty well.
m = hilb(6);
makeHtmlTable(m,[],[],[],[],6);

%%

cmap = jet/2 + 0.5;
makeHtmlTable(m,[],[],[],cmap,2);


%%
m = [1 2; 3 4];
makeHtmlTable(m);

%%
makeHtmlTable(m,[],[],[],summer);

%%
makeHtmlTable(m,[],{'top','bottom'},{'left','right'},winter);

%%

m = magic(11);
makeHtmlTable(m,[],[],[],autumn);

%%
% If both the numeric input is NaN and the string input is non-empty for a
% given cell, then display the string.

m = magic(12);
t = cell(size(m));
m(23) = NaN;
t{23} = 'unknown';
m(79) = NaN;
t{79} = 'missing';
m(120) = NaN;
t{120} = 'no data';
makeHtmlTable(m,t,[],[],winter);


%%
% This displays the frequency with which ASCII characters show up in the
% help entry for the PRINT command.

t = num2cell(char(reshape(48:127,8,10)'));
m = zeros(size(t));
h = help('print');
for i = 1:numel(t)
    m(i) = sum(h==t{i});
end
makeHtmlTable(m,t,[],[],summer);


%% 
% You can even use this to display images. Here's the full image.
load mandrill
image(X);
colormap(map)


%%
% Here's a zoomed sub-image of one eye.
c = 306; 
cinc = 34;
r = 10; 
rinc = 26;
X2 = X((1:rinc)+r, (1:cinc)+c);
image(X2);
colormap(map)

%%

makeHtmlTable(X2,[],[],[],map);
