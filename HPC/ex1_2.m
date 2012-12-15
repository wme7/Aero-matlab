%% Run MATLAB function on multiple cores
j = batch(@prunbirthday, 1, {2e5, 30}, 'configuration', 'hpcfig1',...
    'AttachedFiles',{'prunbirthday.m', 'birthday.m'});

%% Wait for the job to finish
wait(j)   

%% Get results into a cell array
r = fetchOutputs(j)

%% Clean up a batch job's data
delete(j)