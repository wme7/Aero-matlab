%% Create MATLAB pool job
% configuration??????????II. (6)?Profile????
j = createMatlabPoolJob('configuration', 'hpcfig1');
% Directories and files that worker can access
set(j, 'FileDependencies', {'RunDanmIt2.m'})
 
%% Create new task in job
createTask(j, @RunDanmIt2);
% Submit a Job to the Job Queue
submit(j);
% Wait for a job to finish running
waitForState(j);
% Retrieve the Job's Results
results = getAllOutputArguments(j);
results{:}
% Remove job
%destroy(j);