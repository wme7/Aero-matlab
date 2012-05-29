Plot a 3D array using patch
===========================
 
Adam H. Aitkenhead 
adam.aitkenhead@christie.nhs.uk 
The Christie NHS Foundation Trust 
17th August 2010 
 
 
USAGE 
===== 
 
This function enables a 3D array to be displayed using a patch surface mesh.  To plot a 3D logical array, the function is called using the following syntax:
 
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ);
 
Alternatively, a 3D numeric array can be plotted such that the colour of each facet corresponds to the value of each voxel in the array.  The plot is generated using the following command.  (Note that voxels which are not to be displayed should contain a value of NaN.)
 
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,'col');
 
The colormap to be used for the display of a 3D numeric array can be defined as follows:
 
>> cmap = jet(16);
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,cmap,'col');
 
Also, the colorbar lower and upper limits can be defined by the user as follows:
 
>> clim = [1,10];
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,clim,'col');
 
 
INPUT PARAMETERS 
================ 
 
gridINPUT        - 3D array of size (P,Q,R)  - If not using the flag 'col', then gridINPUT should be a logical array.  If using the flag 'col', then gridINPUT should be a numeric array, and all voxels which are not to be displayed should contain a value of NaN.
gridX (optional) - A 1xP array   - List of the X axis coordinates.
gridY (optional) - A 1xQ array   - List of the Y axis coordinates.
gridZ (optional) - A 1xR array   - List of the Z axis coordinates.
cmap  (optional) - A Nx3 array   - The colormap definition.  When plotting using the 'col' flag, cmap must be an Nx3 array, eg jet(32).  When plotting a logical array, cmap must be an RGB triplet, eg [0.5,0.5,0]. 
clim  (optional) - A 2x1 array   - The colormap upper and lower limits.
 
 
ADDITIONAL INPUT FLAGS
======================
 
'col'  (optional) - When this flag is present, the surface is plotted using colours which correspond to the value in each voxel.  In the input array gridINPUT, voxels which are not to be displayed should have a value of NaN.
'barN' (optional) - Display a colorbar on North of plot.
'barE' (optional) - Display a colorbar on East of plot.
'barS' (optional) - Display a colorbar on South of plot.
'barW' (optional) - Display a colorbar on West of plot.

 
 
OUTPUT PARAMETERS 
================ 
 
hpat  (optional)  - Handle to the patch object.
hcbar (optional)  - Handle to the colorbar.
 
 
 
EXAMPLE 
=======
 
For two examples, run the following code:

>> load exampleA.mat
>> figure
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ);

>> load exampleB.mat
>> figure
>> cmap = jet(16);
>> hpat = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,'col',cmap);
 