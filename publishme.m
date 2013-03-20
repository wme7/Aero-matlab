function new_doc=publishme(func,opts,add_fun,rm_fun)
%%
% Publish to not avoid evaluating the published code.
% Mofified by Manuel Diaz, NTU, 2013.03.19.
%
% publishdepfun creates a single published html file with called
% functions attached to the end of the root mfile. Hyperlinks will be 
% within a "Called Functions" in the Table of Contents area of the
% published document that link to the called function.  Functions that are
% called that are under MATLABROOT are not included in the published
% output. 
%
% publishdepfun only goes one level deep.... 
% 
% Example: If fun1 calls fun2, fun2 will be published. If fun2 calls fun3,
%          fun3 will not be published.
%
% Calling: new_doc = publishdepfun(func,opts,add_fun)
% 
% Inputs: 
%
% * func => character string of root function or script ('Function1'). 
%           The ".m" is not needed.
% * opts => is the list of options for the publishing. See publish.m help 
%           for more information.
% * add_fun => is a cell array of additional functions to be published.
% * rm_fun => is a cell array of functions to remove from publish list.
%
%
% Outputs: 
%
% * new_doc => the final html document with attached subfunctions.
%                    Name is the root file name with "_withfuncs" appended.
%
%
%
% Author: Nick Angelini 
%
% Version: 2.1
%
% Date: 11/29/2011
%
%
%

%% 
% Changelog
%
% Version 2.3
%
% # added ability to cell array of functions to remove from publish list
%
%
% Version 2.2
%
% # added ability to remove mex files from output.
%
%
% Version 2.1
%
% # added support for displaying published html document if no output
% variable is specified.
%
%
% Version 2
%
% # added ability to pass additional functions as a cell array to be published.
% # added '-calltree' opition to the depfun call to find the functions
% passed into other functions like ode's and quad. 
%
%
%%
% Known Issues
%
% # Likely dealing with publish options structure isn't honored correctly
% across publishing root and children.
% # Only some basic testing on HTML output, no testing and not expected to
% work for other outputs
% # Can't call publishdepfun on functions that call publish (itself or
% others).  Hit some errors that haven't been explored.
% # Current limitations are if depfun does not return function name it will
% not get published unless the user inputs a secondary list of
% functions. Examples can be found in the help documentation for
% depfun.
% # Currently requires cell blocks to be used, with atleast one block with
% a title other than the overall title of the root file. Need a "Contents"
% section to be created.

%% Error Checking
if nargin<1,error('at least pass in a function'),end
% Default to HTML publishing.
if (nargin < 2)
    opts = 'html';
end
if nargin > 2 && isempty(opts)
    opts = 'html';
end
if nargin < 3
    add_fun = [];
end
if nargin < 4
    rm_fun = [];
end

%% get functions at the top level only
list=depfun(func,'-toponly','-quiet','-calltree');
out=list;

%% Remove hits that are in matlabroot
out(strncmp(matlabroot,list,length(matlabroot)))=[];

%% Remove hits that are probably MEX functions 
idx=strfind(out,mexext); 
for i=1:length(idx), 
    idxs(i) = ~isempty(idx{i}); 
end 
out(idxs)=[];

%% Add user list to functions to be published
out = [out;add_fun];

%% Remove user list
for i=1:length(rm_fun)
    idx=strfind(out,rm_fun{i});
    for i=1:length(idx),
        idxs(i) = ~isempty(idx{i});
    end
    out(idxs)=[];
end

%% Publish main/root/caller function
%out_doc{1}=publish(out{1},opts,);
out_doc{1}=publish(out{1},'format',opts,'evalCode',false);

%% Publish Dependent Functions
options.evalCode=false;
options.outputDir='';

for i=2:length(out)
    out_doc{i}=publish(out{i},options);
end

%% Initialize variables for combining html files
% check is to find id the correct line to add html code, and check2 is to
% find the correct spot in the line to add code for the "Called Functions"
% heading.

check = '</style></head><body>';
check2 = '</div>';
header = sprintf('<h2>Called Functions</h2><ul>');

%% Collect into one HTML

% get name of main html doc
[pathstr, name, ~] = fileparts(out_doc{1});

% create new html with same file name + _withfuncs label
new_doc=[pathstr filesep name '_withfuncs.html'];
fid_new=fopen(new_doc,'w');

% append html docs together
for i=1:length(out_doc)
    % if first/main html
    if i == 1
        % open first html file. not the _withfuncs label
        fid1 = fopen(out_doc{1});

        % loop through line by line of first html doc
        while 1
            x = fgetl(fid1);
            % break out of loop if end of doc is reached
            if x == -1
                break
            end
            % see if current line matches the check variable
            if ~isempty(regexp(x,check,'ONCE'))
                % find location in check line to place subfunction heading
                id = regexp(x,check2);
                % output current line up to char before check2 + header
                counts = fprintf(fid_new,'%s',[x(1:id-1) header]);
                
                % loop through out_doc adding hyperlink code for
                % subfunctions
                for j=2:length(out_doc)
                    % get function name
                    [~, mname, ~] = fileparts(out_doc{j});
                    % store function name for hyperlink flag
                    namestore{j} = mname;
                    % create link code and append to current line
                    newLink = sprintf('<li><a href="#%s">%s </a></li>',mname,mname);
                    count = fprintf(fid_new,'%s',newLink);
                    % when done with hyperlinks append rest of original
                    % line
                    if j == length(out_doc)
                        count = fprintf(fid_new,'</ul>%s\n',x(id:end));
                    end
                end
                % if not line of interest print the line
            else
                counts = fprintf(fid_new,'%s\n',x);
            end
        end
        % close the original driver html
        fclose(fid1);
        % if not first html doc print the hyperlink flag for subfunction then
        % copy the subfunction html to end of _withfuncs.html
    else
        marker = sprintf('<a name="%s"></a>',namestore{i});
        % read in children html publised html documents
        fileStr = fileread(out_doc{i});
        
        % looking for flags used by function grabcode to allow grabcode to 
        % grab only the root function
        matches = regexp(fileStr,'(##### SOURCE BEGIN #####\n.*\n##### SOURCE END #####)','tokens','once');
        newStr = strrep(fileStr,matches{1},'\n');
        count=fprintf(fid_new,'%s\n%s\n',marker,newStr);
    end
end
fclose(fid_new);
if nargout == 0
    web(new_doc)
end