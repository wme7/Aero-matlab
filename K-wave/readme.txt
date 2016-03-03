____________________________________________________________________________

                                  k-Wave

                    A MATLAB toolbox for the time-domain 
                     simulation of acoustic wave fields
____________________________________________________________________________

VERSION INFORMATION
____________________________________________________________________________

Version 1.1.1, Released 9th October 2015
Written by Bradley Treeby, Ben Cox, and Jiri Jaros

Tested using:
   Mac OS X Lion: MATLAB R2011a, R2014a
   Windows 7 64-bit: MATLAB R2008a through to R2014a
   Linux Ubuntu 12.04, 14.04: R2013a, R2014a

Please report bugs and suggestions on http://www.k-wave.org/forum
The toolbox may be downloaded from http://www.k-wave.org/download.php

NOTE: The photoacoustic reconstruction functions kspaceLineRecon and 
kspacePlaneRecon (all toolbox versions) do not work with R2012b.
____________________________________________________________________________

PRODUCT OVERVIEW
____________________________________________________________________________

k-Wave is an open source MATLAB toolbox designed for the time-domain 
simulation of propagating acoustic waves in 1D, 2D, or 3D [1]. The toolbox
has a wide range of functionality, but at its heart is an advanced numerical
model that can account for both linear and nonlinear wave propagation, an 
arbitrary distribution of heterogeneous material parameters, and power law 
acoustic absorption.

The numerical model is based on the solution of three coupled first-order 
partial differential equations which are equivalent to a generalised form 
of the Westervelt equation [2]. The equations are solved using a k-space 
pseudospectral method, where spatial gradients are calculated using a 
Fourier collocation scheme, and temporal gradients are calculated using a
k-space corrected finite-difference scheme. The temporal scheme is exact in
the limit of linear wave propagation in a homogeneous and lossless medium, 
and significantly reduces numerical dispersion in the more general case.

Power law acoustic absorption is accounted for using a linear integro-
differential operator based on the fractional Laplacian [3]. A split-field 
perfectly matched layer (PML) is used to absorb the waves at the edges of 
the computational domain. The main advantage of the numerical model used in 
k-Wave compared to models based on finite-difference time domain (FDTD) 
schemes is that fewer spatial and temporal grid points are needed for 
accurate simulations. This means the models run faster and use less memory. 
A detailed description of the model is given in the k-Wave User Manual and 
the references below.

[1] B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation 
and reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15,
no. 2, p. 021314, 2010. 
[2] B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling 
nonlinear ultrasound propagation in heterogeneous media with power law 
absorption using a k-space pseudospectral method," J. Acoust. Soc. Am., 
vol. 131, no. 6, pp. 4324-4336, 2012.
[3] B. E. Treeby and B. T. Cox, "Modeling power law absorption and 
dispersion for acoustic propagation using the fractional Laplacian," J. 
Acoust. Soc. Am., vol. 127, no. 5, pp. 2741-2748, 2010.
____________________________________________________________________________

INSTALLATION INSTRUCTIONS
____________________________________________________________________________

The k-Wave toolbox is installed by adding the root k-Wave folder to the 
MATLAB path. This can be done using the "Set Path" dialog box which is 
accessed by typing "pathtool" at the MATLAB command line. This dialog box 
can also be accessed using the dropdown menus (File -> Set Path) if using 
MATLAB 2012a and earlier, or the the "Set Path" button on the ribbon bar if
using MATLAB 2012b and later. Once the dialog box is open, the toolbox is 
installed by clicking "Add Folder", selecting the k-Wave toolbox folder, and
clicking "save". The toolbox can be uninstalled in the same fashion. 

For Linux users, using the "Set Path" dialog box requires write access to 
pathdef.m. This file can be found under <...matlabroot...>/toolbox/local. To
find where MATLAB is installed, type "matlabroot" at the MATLAB command line. 
Alternatively, the toolbox can be installed by adding the line 

    addpath('<...pathname...>/k-Wave');
    
to the startup.m file, where <...pathname...> is replaced with the location 
of the toolbox, and the slashes should be in the direction native to your 
operating system. If no startup.m file exists, create one, and save it in 
the MATLAB startup directory. 

After installation, restart MATLAB. You should then be able to see the 
k-Wave help files in the MATLAB help browser. This can be accessed by 
selecting "k-Wave Toolbox" from the contents page. In versions of MATLAB 
prior to 2012b, the help browser is opened by clicking on the blue question
mark icon on the menu bar. In MATLAB 2012b (and later), the documentation is
accessed by selecting "Help" from the ribbon bar, and then clicking on 
"Supplemental Software". Try selecting one of the examples and then clicking
"run the file". 

If you can't see "k-Wave Toolbox" in the contents list of the MATLAB help 
browser, try typing "help k-Wave" at the command prompt to see if the 
toolbox has been installed correctly. If it has and you still can't see the
help files, open "Preferences" and select "Help" and make sure "k-Wave 
Toolbox" or "All Products" is checked. 

After installation, to make the k-Wave documentation searchable from within
the MATLAB help browser, run

    builddocsearchdb(`<...pathname...>/k-Wave/helpfiles');
    
again using the slash direction native to your operating system. Note, the
created database file will only work with the version of MATLAB used to 
create it.

If using the C++ or CUDA versions of kspaceFirstOrder3D, the appropriate 
binaries (and library files if using Windows) should also be downloaded 
from http://www.k-wave.org/download.php and placed in the root "binaries"
folder of the toolbox.
____________________________________________________________________________

RELEASE NOTES V1.1.1
____________________________________________________________________________

Bug Fixes:
- bug fix in pstdElastic2D and pstdElastic3D when defining material 
  coefficients on staggered grid (generated incorrect results)
____________________________________________________________________________

RELEASE NOTES V1.1
____________________________________________________________________________

New Features and Changes:
- simulations in elastic media are now supported using the pstdElastic2D and
  pstdElastic3D functions
- 3D simulations can be run on an NVIDIA graphics processing unit (GPU) 
  using a native C++/CUDA code
- rectangular binary sensor masks can be defined by assigning the grid 
  coordinates of two opposing corners to sensor.mask
- sensor.record now supports the additional input options 'p_min', 
  'p_max_all', 'p_min_all', 'u_min', 'u_max_all', 'u_min_all', and 
  'u_non_staggered'
- simulation functions have been restructured to reduce code duplication
- the accuracy of the acoustic intensity output for frequencies close to 
  Nyquist limit has been improved
- the C++ code now supports command line inputs for checkpoint-restart and 
  importing FFTW wisdom
- the C++ code now exploits AVX instructions in addition to SSE
- attenuation compensation for photoacoustic tomography can be performed 
  using time-variant filtering via attenComp
- vessel filtering for photoacoustic tomography can be performed using 
  vesselFilter
- kspaceFirstOrder3DC now supports optional inputs to specify the name and 
  location of the binary and data files
- makeSphere and makeSphericalSection can be set to return logical matrices 
  to reduce memory usage
- resize now supports 1D inputs
- checkFactors now also reports prime numbers

Bug Fixes:
- bug fix in using a Cartesian sensor mask with nearest neighbour 
  interpolation and sensor.record set to 'p_min' (generated error)
- bug fix in using sensor.directivity_angle with the optional input 'DataCast'
  set to 'gpuArray-single' (generated error)
- bug fix in kspaceFirstOrder3DC when performing time reversal image 
  reconstruction using sensor.time_reversal_boundary_data with 'PMLInside' set
  to false (PML was not removed)
- bug fix in kspaceFirstOrder3DC when calling script on a different drive 
  using Windows (generated error)
- bug fix in kspaceSecondOrder when setting absorption values with 
  medium.alpha_power to values less than 1 (returned NaN)
- bug fix in kspaceSecondOrder when using sensor.record = {'p_final'} 
  (generated error, grid expansion not removed)
- bug fix in makeSphere for grid sizes with uneven dimensions (generated 
  error)
- bug fix in makeSphere to centre sphere on the kgrid origin (previously 
  offset by 1)
- bug fix in writeMatrix for source variables larger than 1 MB (incorrect 
  chunk size generated error)

New Functions:
- attenComp
- checkStability
- fitPowerLawParams
- getBLI
- getSpacedPoints
- h5compare
- kspaceFirstOrder3DG
- maxND
- minND
- overlayPlot
- pstdElastic2D
- pstdElastic3D
- revolve2D
- timeShift
- vesselFilter
- writeAttributes
- writeFlags
- writeGrid

New Examples:
- Defining A Sensor Mask By Opposing Corners
- Attenuation Compensation Using Time Variant Filtering
- Running C++ Simulations
- Saving Input Files in Parts
- Explosive Source In A Layered Medium
- Plane Wave Absorption
- Shear Waves And Critical Angle Reflection
- Simulations In Three Dimensions
____________________________________________________________________________

LICENSE
____________________________________________________________________________

k-Wave (c) 2009-2014 Bradley Treeby and Ben Cox

The k-Wave toolbox is distributed by the copyright owners under the terms of
the GNU Lesser General Public License (LGPL) which is a set of additional 
permissions added to the GNU General Public License (GPL). The full text of 
both licenses is included with the toolbox in the folder 'license'.

The licence places copyleft restrictions on the k-Wave toolbox. Essentially,
anyone can use the software for any purpose (commercial or non-commercial), 
the source code for the toolbox is freely available, and anyone can 
redistribute the software (in its original form or modified) as long as the
distributed product comes with the full source code and is also licensed 
under the LGPL. You can make private modified versions of the toolbox 
without any obligation to divulge the modifications so long as the modified
software is not distributed to anyone else. The copyleft restrictions only 
apply directly to the toolbox, but not to other (non-derivative) software 
that simply links to or uses the toolbox. 

k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details (http://www.gnu.org/licenses/lgpl.html). 

If you find the toolbox useful for your academic work, please consider 
citing:

B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and 
reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15, no. 
2, p. 021314, 2010.

and/or

B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear 
ultrasound propagation in heterogeneous media with power law absorption 
using a k-space pseudospectral method," J. Acoust. Soc. Am., vol. 131, 
no. 6, pp. 4324-4336, 2012.

along with any other relevant publications. The first paper gives an 
overview of the toolbox with applications in photoacoustics, and the second 
describes the nonlinear ultrasound model and the C++ code. 
____________________________________________________________________________