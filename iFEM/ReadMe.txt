INSTALLATION

Add the path to iFEM into the path library of MATLAB. There are two ways: 
1. Graphical interface. Click
	File -> Set Path -> Add with Subfolders
   and chose the directory where the package iFEM is stored.
2. Command window. Go to the directory of iFEM and run 
>> setpath
Note: if you want to include the path 
------------------------------------------------------------------
HELP

1. HELP FUN displays a description of and syntax for the function FUN. For example,
>> help mg
will show basic usage for mg function in the plain text.  

2. IFEM FUNdoc show detailed description. For example,
>> ifem mgdoc
will explain the mg function step by step in html format.

------------------------------------------------------------------
QUICK START

1. Type 
>>ifem introduction 
to get an introduction on ifem.

2. Go through examples in \example directory.

------------------------------------------------------------------
FEED BACK

If you like it, please send me an email lyc102@gmail.com. If you  feel it is helpful for your research, please acknowledge your use by citing:
 
L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, University of California at Irvine, 2009.

------------------------------------------------------------------
ACKNOWLEDGEMENT

The author would like to thank Professor Michael Holst in University of California at San Diego and Professor Ludmil Zikatanov in Pennsylvania State University for many insightful discussion, and also Dr. Chensong Zhang in Pennsylvania State University for the effort in the development of AFEM@matlab, an early version of iFEM.

The author also thanks students Ming Wang and Huayi Wei for their contribution to iFEM in one way or another. Detailed credits can be found in 
the M-lint of m files.


Long Chen

--------------------------------------------------------------------
Associate Professor                   Fax: 949-824-7993
Department of Mathematics             Phone: 949-824-6595
University of California at Irvine    Email: chenlong@math.uci.edu
Irvine, CA, 92697                     http://math.uci.edu/~chenlong/
--------------------------------------------------------------------

