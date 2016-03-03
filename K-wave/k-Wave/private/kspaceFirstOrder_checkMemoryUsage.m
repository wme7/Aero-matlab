% DESCRIPTION:
%       subscript to extract and display GPU and CPU memory usage
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 16th July 2013
%       last update - 16th July 2013
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% display current matlab memory usage
if nargout == 2 && strncmp(computer, 'PCWIN', 5)
    [mem_usage.user, mem_usage.sys] = memory;
    disp(['  memory used: ' num2str(mem_usage.user.MemUsedMATLAB./1024^3) ' GB (of ' num2str(mem_usage.sys.PhysicalMemory.Total./1024^3) ' GB)']); 
end        

% gpu memory counter for GPUmat toolbox
if strncmp(data_cast, 'kWaveGPU', 8)
    current_gpu_mem = GPUmem;
    disp(['  GPU memory used: ' num2str((total_gpu_mem - current_gpu_mem)/2^30) ' GB (of ' num2str(total_gpu_mem/2^30) ' GB)']);
    mem_usage.gpu.total = total_gpu_mem;
    mem_usage.gpu.used = total_gpu_mem - current_gpu_mem;            
end

% gpu memory counter for Accelereyes toolbox
if strcmp(data_cast, 'gsingle') || strcmp(data_cast, 'gdouble')
    gpu_info = ginfo(true);
    disp(['  GPU memory used: ' num2str((gpu_info.gpu_total - gpu_info.gpu_free)/2^30) ' GB (of ' num2str(gpu_info.gpu_total/2^30) ' GB)']);
    mem_usage.gpu.total = gpu_info.gpu_total;
    mem_usage.gpu.used = gpu_info.gpu_total - gpu_info.gpu_free;              
end   

% gpu memory counter for Parallel Computing toolbox
if strcmp(data_cast, 'gpuArray')
    gpu_info = gpuDevice;
    disp(['  GPU memory used: ' num2str((gpu_info.TotalMemory - gpu_info.FreeMemory)/2^30) ' GB (of ' num2str(gpu_info.TotalMemory/2^30) ' GB)']);
    mem_usage.gpu.total = gpu_info.TotalMemory;
    mem_usage.gpu.used = gpu_info.TotalMemory - gpu_info.FreeMemory;            
end