classdef PerformanceData
    %PERFORMANCEDATA  a class to store GPUBench performance data
    %
    %   p = gpuBench measures some performance data for your currently
    %   selected GPU. Do not use this class directly - always use gpuBench
    %   to create the data.
    %
    %   See also: gpuBench
    
    %   Copyright 2011 The MathWorks, Inc.
    
    properties
        IsSelected
        IsHostData
    end
    
    properties (SetAccess=private)
        MATLABRelease
        CPUInfo
        GPUInfo
        Results
        Timestamp
    end
    
    methods
        function obj = PerformanceData( release, cpuinfo, gpuinfo, timestamp )
            % Construct a new performance data object
            obj.MATLABRelease = release;
            obj.CPUInfo = cpuinfo;
            obj.GPUInfo = gpuinfo;
            obj.Results = struct( ...
                'FunctionName', {}, ...
                'DataType', {}, ...
                'Sizes', {}, ...
                'NumOps', {}, ...
                'Times', {} );
            obj.IsSelected = false;
            obj.IsHostData = false;
            obj.Timestamp = timestamp;
        end % constructor
        
        function obj = addResult( obj, fcnName, datatype, sizes, numops, times )
            N = numel( obj.Results );
            obj.Results(N+1,1).FunctionName = fcnName;
            obj.Results(N+1,1).DataType = datatype;
            obj.Results(N+1,1).Sizes = sizes;
            obj.Results(N+1,1).NumOps = numops;
            obj.Results(N+1,1).Times = times;
        end % addResult
        
        function out = hasResult( obj, fcnname, datatype )
            if nargin<3
                % Name may be 'fcn (type)'
                [fcnname,datatype] = iSplitName( fcnname );
            end
            nameMatches = ismember( {obj.Results.FunctionName}, fcnname );
            typeMatches = ismember( {obj.Results.DataType}, datatype );
            out = any( nameMatches & typeMatches );
        end % hasResult
        
        function out = getResult( obj, fcnname, datatype )
            if nargin<3
                % Name may be 'fcn (type)'
                [fcnname,datatype] = iSplitName( fcnname );
            end
            nameMatches = ismember( {obj.Results.FunctionName}, fcnname );
            typeMatches = ismember( {obj.Results.DataType}, datatype );
            idx = find( nameMatches & typeMatches, 1, 'first' );
            if isempty( idx )
                error( 'GPUBench:PerformanceData:NoSuchData', 'No results were found for %s (%s).', ...
                    fcnname, datatype );
            end
            out = obj.Results(idx);
        end % getResult
        
        function name = getDeviceName( obj )
            if obj.IsHostData
                name = 'Host PC';
            else
                name = obj.GPUInfo.Name;
            end
        end
    end
    
end

function [fcnname,datatype] = iSplitName( longname )
% Split a long name 'fcn (datatype)' into its component name and type
out = regexp( longname, '(?<fcn>\w+)\s+\((?<type>\w+)\)', 'names' );
fcnname = out.fcn;
datatype = out.type;
end
