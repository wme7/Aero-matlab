%% Arrayfun Example
clear all; close all; clc;

%% Example
% Run a function using the GPU by asking matlab to evaluate this function
% with data created inside our Gpu memory (gpuarray type data)
D = gpuDevice(1); % because, I have 2 GPUs in my PC.

% Create gpuarrays data,
s1 = gpuArray(rand(400));
s2 = gpuArray(rand(400));
s3 = gpuArray(rand(400));

% Evaluate aGpuFunction with gpuarray arguments,
[o1,o2] = arrayfun(@aGpuFunction,s1,s2,s3);

% We get a gpurray solution, therefore we need to gather it from gpuDevice
d = gather(o2);

% See variables
whos

% Clear GPU Memory
reset(D);
