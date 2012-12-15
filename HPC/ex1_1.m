%% Create MATLAB pool job
j = createMatlabPoolJob('configuration', 'hpcfig1');

% Directories and files that worker can access
set(j, 'FileDependencies', {'prunbirthday.m', 'birthday.m'})

%% Create new task in job
createTask(j, @prunbirthday, 1, {2e5, 30});

%% Submit a Job to the Job Queue
submit(j);

%%  Wait for a job to finish running
waitForState(j);

%% Retrieve the Job's Results
results = getAllOutputArguments(j);
results{:}

%% Remove job
destroy(j);

