function match = birthday(groupsize)
% BIRTHDAY Simulates a single trial of the Birthday Paradox
%    MATCH = BIRTHDAY(GROUPSIZE) creates a randomly-selected birthday for
%    every member of a group of size GROUPSIZE and tests whether any of 
%    the selected birthdays match.  MATCH is 1 if two or more members of 
%    the group share the same birthday and 0 otherwise.
%
%    Example:
%    >> match = birthday(30)


% Match is zero until a birthday match is found
match = 0;

% Initialize list of taken birthdays
bdaylist = zeros(1, groupsize);

for person = 1:groupsize

    % Randomly select a birthdate for the individual (ignore leap years)
    birthdate = ceil(365*rand());

    % Check if someone else in the group shares the same birthday
    if (any(birthdate == bdaylist))

        % A match is found, return from the function
        match = 1;
        return;
        
    end

    % Add the birthdate to the list for the group
    bdaylist(person) = birthdate;

end
    