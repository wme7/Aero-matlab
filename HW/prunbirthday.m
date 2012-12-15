function prob = prunbirthday(numtrials, groupsize)
% PRUNBIRTHDAY Runs a parallel version of the Birthday Paradox simulation
% PROB = PRUNBIRTHDAY(NUMTRIALS, GROUPSIZE) Calls the birthday code
% NUMTRIALS times to see if any birthdays match in a group of size GROUPSIZE.  
% The return value is the probability that a match will be found.

% Preallocate some memory for the matches
matches = zeros(1, numtrials);

parfor trial = 1:numtrials
    % Run a simulation for a group
    matches(trial) = birthday(groupsize);
% b(trial) = matches(trial-1); dependent case
end

% Probability is the sum of matches divided by number of trials
prob = sum(matches)/numtrials;
