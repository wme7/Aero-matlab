%% Trying 2 different solvers on the same ODE system
% Looking at Spring Damper System and comparing the results
% from the ODE45, ODE23, ODE15s solvers with the analytical solution.
%
% Copyright 2009-2011 The MathWorks, Inc.

%% Setting Initial Conditions
m = 5;  % mass
b = linspace(0.1, 5, 40);  % damping value
k = linspace(1.5, 5, 40);  % stiffness value
totalTime = 25; % total time of simulation

%% Compute Serially

disp('Computing serially...');
tic;

p1 = springDampSolver15s(m, k, b, totalTime);   % ODE15s
p2 = springDampSolver23(m, k, b, totalTime);    % ODE23
p3 = springDampSolver45(m, k, b, totalTime);    % ODE45
p4 = springDampSolverAna(m, k, b, totalTime);   % Analytical

% Compare results
compareResults(b, k, p1, p2, p3, p4)

t1 = toc;
fprintf('\nElapsed time is %0.2f seconds.\n', t1);

%% Compute Using Jobs and Tasks

disp('Computing using Jobs and Tasks...');
tic;

% Find resource/scheduler
sched = findResource();

% Create job
job = createJob(sched);
set(job, 'FileDependencies', ...
   {'springDampSolver15s.m', 'springDampSolver23.m', ...
   'springDampSolver45.m', 'springDampSolverAna.m', 'odesystem.m'})

% Create tasks
task1 = createTask(job, @springDampSolver15s, 1, {m, k, b, totalTime});
task2 = createTask(job, @springDampSolver23 , 1, {m, k, b, totalTime});
task3 = createTask(job, @springDampSolver45 , 1, {m, k, b, totalTime});
task4 = createTask(job, @springDampSolverAna, 1, {m, k, b, totalTime});

% Submit job
submit(job)
wait(job) % optional, only load when job is finished

% Retrieve results
results = getAllOutputArguments(job);

% Compare results
compareResults(b, k, results{1,1}, results{2,1}, ...
   results{3,1}, results{4,1});

t2 = toc;
fprintf('\nElapsed time is %0.2f seconds.\n', t2);

% Destroy job when finished
destroy(job)

%% Speed Up
% When you manually create tasks using Jobs and Tasks, you need to be
% aware that the way you split up the tasks will affect the performance. In
% this example, the various solvers take different amount of time to solve
% the same problem, so you may not get optimal utilization of your
% resources. Having similar-sized computations or large number of tasks may
% help even out the work load.

fprintf('\n\nSpeed up (time serial / time parallel): %0.2f\n', t1/t2);