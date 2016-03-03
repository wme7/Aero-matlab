% k-Wave Toolbox
% Version 1.1 03-Oct-2014
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See k-Wave Toolbox in help menu for description and examples.
% Type "help <command-name>" for documentation on individual commands.
% -----------------------------------------------------------------
%
% Wave Propagation
%   kspaceFirstOrder1D      - 1D time-domain simulation of wave propagation
%   kspaceFirstOrder2D      - 2D time-domain simulation of wave propagation
%   kspaceFirstOrder3D      - 3D time-domain simulation of wave propagation
%   kspaceFirstOrder3DC     - 3D time-domain simulation of wave propagation using C++ code
%   kspaceFirstOrder3DG     - 3D time-domain simulation of wave propagation on a GPU using C++ CUDA code
%   kspaceSecondOrder       - Fast time-domain simulation of wave propagation for homogeneous media
%   mendousse               - Compute Mendousse's solution for nonlinear wave propagation in viscous media
%   pstdElastic2D           - 2D time-domain simulation of elastic wave propagation
%   pstdElastic3D           - 3D time-domain simulation of elastic wave propagation
% 
% Image Reconstruction
%   kspaceLineRecon         - 2D linear FFT reconstruction
%   kspacePlaneRecon        - 3D planar FFT reconstruction
%
%   See also kspaceFirstOrder1D, kspaceFirstOrder2D, and kspaceFirstOrder3D for time-reversal image reconstruction
% 
% Geometry Creation
%   makeBall                - Create a binary map of filled ball within a 3D grid
%   makeCartCircle          - Create a 2D Cartesian circle or arc
%   makeCartSphere          - Create a 3D Cartesian sphere
%   makeCircle              - Create a binary map of a circle within a 2D grid
%   makeDisc                - Create a binary map of a filled disc within a 2D grid
%   makeLine                - Create a binary map of a straight line within a 2D grid
%   makeSphere              - Create a binary map of a sphere within a 3D grid
%   makeSphericalSection	- Create a binary map of a sphere segment within a 3D grid
%
% Acoustic Absorption
%   attenComp               - Attenuation compensation using time-variant filtering
%   attenuationWater        - Calculate ultrasound attenuation in distilled water
%   db2neper                - Convert decibels to nepers
%   fitPowerLawParams       - Fit power law absorption parameters for highly absorbing media
%   neper2db                - Convert nepers to decibels
%   powerLawKramersKronig	- Calculate dispersion for power law absorption
%
% Grid and Matrix Utilities
%   cart2grid               - Interpolate a set of Cartesian points onto a binary grid
%   expandMatrix            - Enlarge a matrix by extending the edge values
%   findClosest             - Return the closest value in a matrix
%   grid2cart               - Return the Cartesian coordinates of the non-zero points of a binary grid
%   interpCartData          - Interpolate data from a Cartesian to a binary sensor mask
%   interpftn               - Resample data using Fourier interpolation
%   loadImage               - Load an image file
%   makeGrid                - Create k-Wave grid structure
%   maxND                   - Return the value and indices of the largest value in an N-D array
%   minND                   - Return the value and indices of the smallest value in an N-D array
%   numDim                  - Return the number of matrix dimensions
%   resize                  - Resize a matrix
%   reorderSensorData       - Reorder sensor data from kspaceFirstOrder2D based on angle
%   revolve2D               - Form 3D matrix from revolution of 2D matrix
%   timeShift               - Shift time series to and from staggered temporal grid
%   unmaskSensorData        - Reorder data recorded using a binary sensor mask
%
% Filtering and Spectral Utilities
%   applyFilter             - Filter input with low, high, or band pass filter
%   envelopeDetection       - Extract signal envelope using the Hilbert Transform
%   filterTimeSeries        - Filter signal using the Kaiser windowing method
%   gaussianFilter          - Filter signals using a frequency domain Gaussian filter
%   getAlphaFilter          - Create filter for medium.alpha_filter
%   getBLI                  - Compute underlying Fourier band-limited interpolant (BLI)
%   getFDMatrix             - Create a matrix of finite-difference coefficients
%   getWin                  - Return a frequency domain windowing function
%   gradientFD              - Calculate the gradient using a finite-difference method
%   gradientSpect           - Calculate the gradient using a Fourier spectral method
%   sharpness               - Calculate image sharpness metric
%   spect                   - Compute the single sided amplitude and phase spectrums
%   smooth                  - Smooth a matrix
%   vesselFilter            - Frangi's 3D vessel filter
%
% Display and Visualisation
%   beamPlot                - Plot volumetric data using intersecting planes
%   flyThrough              - Display a three-dimensional matrix slice by slice
%   getColorMap             - Return default k-Wave color map
%   overlayPlot             - Overlay two images
%   saveTiffStack           - Save volume data as a tiff stack
%   scaleFig                - Resize current figure window
%   stackedPlot             - Stacked linear plot
%   voxelPlot               - 3D plot of voxels in a binary matrix
%
% Ultrasound Utilities
%   envelopeDetection       - Extract signal envelope using the Hilbert Transform
%   focus                   - Create input signal based on source mask and focus position
%   gaussianFilter          - Filter signals using a frequency domain Gaussian filter
%   hounsfield2density      - Convert Hounsfield units to density
%   logCompression          - Log compress an input signal
%   makeTransducer          - Create k-Wave ultrasound transducer
%   scanConversion          - Convert scan-lines in polar coordinates to a B-mode ultrasound image
%   toneBurst               - Create an enveloped single frequency tone burst
%
% HDF5 Utilities
%   h5compare               - Compare the contents of two HDF5 files
%   writeAttributes         - Write attributes to a k-Wave HDF5 file
%   writeFlags              - Write input flags to a k-Wave HDF5 file
%   writeGrid               - Write grid and PML properties to a k-Wave HDF5 file
%   writeMatrix             - Write MATLAB matrix to a k-Wave HDF5 file
%
% Other Utilities
%   addNoise                - Add Gaussian noise to a signal for a given SNR
%   benchmark               - Run performance benchmark
%   checkFactors            - Return the maximum prime factor for a range of numbers
%   checkStability          - Return maximum stable timestep for kspace models
%   fwhm                    - Compute the full width at half maximum
%   gaussian                - Create a Gaussian distribution
%   getDateString           - Create a string of the current date and time
%   getkWavePath            - Return pathname to the k-Wave Toolbox
%   getSpacedPoints         - Create vector of log or linear spaced points
%   makeTime                - Create an evenly spaced array of time points
%   scaleSI                 - Scale a number to nearest SI unit prefix
%   scaleTime               - Convert seconds to hours, minutes, and seconds
%   speedSoundWater         - Calculate the speed of sound in water with temperature
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
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