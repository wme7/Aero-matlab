function git(varargin)
% A thin MATLAB wrapper for Git.
% 
%   Short instructions:
%       Use this exactly as you would use the OS command-line verison of Git.
% 
%   Long instructions are:
%       This is not meant to be a comprehensive guide to the near-omnipotent 
%       Git SCM:
%           http://git-scm.com/documentation
% 
%       Common MATLAB workflow: 
% 
%       % Creates initial repository tracking all files under some root
%       % folder
%       >> cd ~/
%       >> git init
%
%       % Shows changes made to all files in repo (none so far)
%       >> git status
%
%       % Create a new file and add some code
%       >> edit foo.m
%
%       % Check repo status, after new file created
%       >> git status
%
%       % Stage/unstage files for commit
%       >> git add foo.m          % Add file to repo or to stage
%       >> git reset HEAD .       % To unstage your files from current commit area
%
%       % Commit your changes to a new branch, with comments
%       >> git commit -m 'Created new file, foo.m'
% 
%       % Other useful commands (replace ellipses with appropriate args)
%       >> git checkout ...       % To restore files to last commit
%       >> git branch ...         % To create or move to another branch
%       >> git diff ...           % See line-by-line changes 
%
%   Useful resources:
%       1. GitX: A visual interface for Git on the OS X client
%       2. Github.com: Remote hosting for Git repos
%       3. Git on Wikipedia: Further reading 
% 
% v0.1,     27 October 2010 -- MR: Initial support for OS X & Linux,
%                               untested on PCs, but expected to work
% 
% v0.2,     11 March 2011   -- TH: Support for PCs
% 
% v0.3,     12 March 2011   -- MR: Fixed man pages hang bug using redirection
% 
% Contributors: (MR) Manu Raghavan
%               (TH) Timothy Hansell


% Test to see if git is installed
[status,~] = system('git status');
% if git is in the path this will return a status of 0 or 128
% depending on whether we are sitting in a repository or not
% it will return a 1 only if the command is not found

    if (status==1)
        % If GIT Is NOT installed, then this should end the function.
        fprintf('git is not installed\n%s\n',...
               'Download it at http://git-scm.com/download');
    else
        % Otherwise we can call the real git with the arguments
        arguments = parse(varargin{:});  
        if ispc
          prog = 'type'
        else
          prog = 'cat'
        end
        [~,result] = system(['git ',arguments,' | ',prog]);
        disp(result)
    end
end

function space_delimited_list = parse(varargin)
    space_delimited_list = cell2mat(...
                cellfun(@(s)([s,' ']),varargin,'UniformOutput',false));
end
