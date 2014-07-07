function [s] = mat2tecplot(tdata,tecfile)
%
% function s =mat2tecplot(tdata,tecfile)
% 
% Program to create tecplot binary data file (.plt).
%
% Format is described here:
%     http://download.tecplot.com/360/dataformat.pdf
%     on page 149. It is also excerpted into text file in
%     ./doc/binaryformat_2012.txt
%
% This program benefitted from the following programs
%   http://tecplottalk.com/viewtopic.php?t=879&sid=2775b1aff4a1dc523ee890b85ab08f40
%
% Here are general rules for binary data in tecplot
%
% http://download.tecplot.com/tecio/360/binaryformat.txt
%
%  1) All character strings in Fortran or matlab must be terminated with
%      null character. In Fortran this is done by "MyString"//char(0)
%      In matlab, this is done by (after writing the string) :
%           dummy_int32 =0;
%           fwrite(fid_out,dummy_int32,'int32');
%
%  2) All boolean Flags  use value 1 for true and 0 for false
%
% A tecplot structure tdata (first argument) is a structure and it
% contains the following information :
%
%     tdata.title    ---title of the plot, string
%                       by default, the title name will be 'tecplot data'
%
%     tdata.Nvar     ---integer, storing number of variables
%
%     tdata.varnames ---cell array of variable names, if array size less than
%                       Nvar, then the rest of them will be named
%                       as var1,var2,var3,... depending on
%                       how many of them have names given. If
%                       length of varnames is larger than Nvar, those 
%                       beyond Nvar will be ignored. varnames{1:Nvar}
%
%     tdata.vformat  ---integer vector array storing dataformat for
%                       each variable (1:Nvar). 1- for float
%                                               2- for double
%                                               3- for longInt
%                                               4- for shortInt
%                                               5- for Byte
%                                               6- for Bit
%                       Default is 1
%
%     tdata.lines    ---an OPTIONAL structure array [Nlines] of line values Y, Z, v 
%                     defined on 1D X ( or with Y, Z for 3D line)
%                     coordinates with I-ordered format. Each line is 
%                     treated as a zone. Lines contain the 
%                     following metadata information (iline=1:Nlines):
% 
%                     lines(iline).varloc --- location of variables, 0 for nodal
%                                            1 for center. 
%                     lines.zonename ---zone names of each lines, if not given
%                                      it will be named as line1,line2,... 
%                     lines.x        ---1D x values of the line (must have)
%                     lines.y        ---1D y values of the line (optional)
%                                       default to zero if not given
%                     lines.z        ---1D z values of the line (optional)
%                                       default to zero, if not given
%                     lines.v        ---2D variables defined on the line
%                                       (optional), default to zero, if not
%                                       given and Nvar>3.  Each variable
%                                       uses one row. E.g. lines.v(2,:) 
%                                       is the 2nd variable in v. Length of
%                                       the row (number of columns) must
%                                       be same as x, y, z. Nvar=3+number of
%                                       rows in v. 
%                     lines.varmin --- minimum value of the line, automatically
%                                      calculated if not given
%                     lines.varmax --- maximum value of the line, automatically 
%                                      calculated if not given
%                     lines.auxname --- auxiliary data name for each line, set to 
%                                       none if not given or empty
%                                       (cell array of strings or single 
%                                        string)
%                     lines.auxval  --- auxiliary data value for each line, set to
%                                       none if not given or empty. Only
%                                       useful if auxname is given.
%                                       (cell array of strings or single 
%                                        string)
%                     lines.datapacking ---is it block or not, by default
%                                          it is point, but block will save
%                                          space. This is deprecated. 
%                                          By default it is set to zero
%                                          (block)
%                     lines.strandID    ---strandID for each of the line
%                                          zone. Each line is treated as a
%                                          zone, and you can assign
%                                          strandID to it. Multiple
%                                          zones can have the same strandID
%                                          such that they grouped together
%                                          and identifiable by the strandID
%                                          If not provided,default is -2
%                                          so that tecplot will assign one 
%                                          automatically. Values can be:
%                                          -2 = strand ID to be assigned by tecplot.
%                                          -1 = static strand ID
%                                          0 <= value < 32700 is valid strand ID
%                     lines.solutiontime ---time associated with the line
%                                           multiple zones can have the
%                                           same solutiontime, hence 
%                                           useful for grouping zones
%                                           together for animation.
%                                           (Optional), default = 0
% 
%     tdata.surfaces ---an array of surface values v and Z defined on 
%                       2D XY coordinates with IJ-ordered format, or v and 
%                       X defined on 2D (YZ) coordinates (JK-ordered 
%                       format) or v and Y defined on 2D (XZ) coordinates 
%                       (IK-ordered format). Each surface is treated
%                       as a zone. The surfaces contain the following 
%                       information to define the zones.
%
%                  surfaces.varloc  --- location of variables for each 
%                                       surface variable, 0-nodal, 1-center
%                                       Note that coordinate variables
%                                       must be nodal. Which variables
%                                       are coordinate variables is
%                                       determined from surfaces.order
%
%                  surfaces.x       --- 2D array giving
%                                       x values of the surface
%                                       (must have if order=1,2)
%                                       (optional if order=1, default to zeros)%
%
%                  surfaces.y       --- 2D array giving 
%                                       y values of the surface
%                                       (must have if order=1,3)
%                                       (optional if order=2, default to zeros)
%
%                  surfaces.z       --- 2D array giving 
%                                       z values of the surface
%                                       (must have if order=2,3)
%                                       (optional if order=3, default to zeros)
%
%                  surfaces.v       ----2D variables defined on the
%                                       surface (optional). Default to
%                                       zeros if Nvar>3
%
%                  surfaces.order   --- determine order (default is 3)
%                                       3- IJ order, surface is z=f(x,y),
%                                          x and y must be nodal. 
%                                          If z not available, default
%                                          to zero
%                                       2- IK order, surface is
%                                          y=f(x,z), x and z must be
%                                          nodal.
%                                          If y not available, default
%                                          to zero
%                                       1- JK order, surface is
%                                          x=f(y,z), y and z must be
%                                          nodal. 
%                                          if x not available, default
%                                          to zero
%
%                  surfaces.zonename --- zone name of each surface, if not
%                                       given, then it will be named as 
%                                       surface1,surface2,...
%                  surfaces.auxname --- auxiliary data name for each surface
%                                        (cell array of strings or single 
%                                        string)
%                  surfaces.auxval  --- auxiliary data value for each auxname
%                                        (cell array of strings or single 
%                                        string)
%                  surfaces.datapacking --is it block or point (deprecated
%                                        and set to zero by default)
%                  surfaces.strandID --- strandID for each surface as 
%                                        a zone. 
%                                        If not provided,default is -2
%                                        so that tecplot will assign one
%                                        automatically. Values can be:
%                                        -2 = strand ID to be assigned by tecplot.
%                                        -1 = static strand ID
%                                        0 <= value < 32700 is valid strand ID
%                  surfaces.solutiontime --- time associated with the 
%                                            surface. (Optional) default 0
%                                           
%     tdata.cubes    ---an array of 3D values v defined on 3D XYZ coordinates
%                      with IJK-ordered format. Each volume dataset is
%                      treated as a zone. cubes contain the following
%                      information to define the zones.
% 
%                      cubes.varloc  --- location of variable for each cube data
%                                        0-nodal, 1-center (default is 0)
%                      cubes.x       --- 3D array giving x coordinate
%                                        values of the surface. (must have)
%                                        size [Imax,Jmax,Kmax]
%                      cubes.y       --- 3D array giving y coordinate values
%                                        of the surface (must have)
%                                        size [Imax,Jmax,Kmax]
%                      cubes.z       --- 3D array giving z coordinate values
%                                        of the surface (must have)
%                                        size [Imax,Jmax,Kmax]
%                      cubes.v       --- 4D array giving value of the variable
%                                        (optional). Default to zeros if
%                                        not given and Nvar>3. First
%                                        dimension is variables, 2nd, 3rd,
%                                        and fourth dimensions are (x,y,z)
%                                        for varilable iv-3, where
%                                        iv=1:Nvar. For nodal data, v is
%                                        of size [Nvar-3,Imax,Jmax,Kmax]; for 
%                                        cell centered data (varloc==1)
%                                        v is of size
%                                        [Nvar-3,Imax-1,Jmax-1,Kmax-1]
%                      cubes.zonename --- zone name of each cube data v, if not
%                                         given, then it will be named as
%                                         cube1,cube2,cube3,... (Optional)
%                      cubes.auxname --- auxiliary data names for each cube
%                                        (cell array of strings or single 
%                                        string) (Optional)
%                      cubes.auxval  --- auxiliary data value for each auxname
%                                        (cell array of strings or single 
%                                        string) (Optional)
%                      cubes.datapacking ---packing method of the cube
%                                           variables (0) for block, 1 for 
%                                           point. Note: for point, varloc
%                                           can't be 1-center. (Deprecated
%                                           and set to zero by default)
%                      cubes.strandID --- standID associated with the cube.
%                                         (Optional).
%                                          If not provided,default is -2
%                                          so that tecplot will assign one
%                                          automatically. Values can be:
%                                          -2 = strand ID to be assigned by tecplot.
%                                          -1 = static strand ID
%                                          0 <= value < 32700 is valid strand ID
%                      cubes.solutiontime ---solution time associated with
%                                            the cube zone (Optional).
%                                            Default is 0
%                      
%     tdata.FElines  ---an array of finite element 1D data, which is                      
%                       a set of line segments defining a 2D or 3D line.
%                       Each finite element line (segments) is treated as a
%                       zone. FElines contain the following information to
%                       define a zone:
% 
%                       FElines.e2n     ---element (made of 2 points) node 
%                                          connectivity list, for each element
%                                          it consists two node numbers, 
%                                          array of size [NE,2] defining
%                                          line element to node connectivity
%                                          where NE is number of elements
%                       FElines.v       ---variable value of each line defined
%                                          on line element (nodal or center
%                                          of line element). 2D dimensional
%                                          1st dimension is variable 
%                                          2nd dimension is nodal or cell
%
%                       FElines.x      ---x coordinate  (1D array)
%                       FElines.y      ---y coordinate  (1D array)
%                       FElines.z      ---z coordinate  (1D array)
%
%                       FElines.order  ---order of each line (default is 4)
%                                         0 : 1D line defined on x
%                                                y=f(x), z=f(x),v=f(x)
%                                         3 : 2D line defined on (x,y)
%                                                z=f(x,y),v=f(x,y)
%                                         2 : 2D line defined on (x,z)
%                                                y=f(x,z),v=f(x,z)
%                                         1:  2D line defined on (y,z)
%                                                x=f(y,z),v=f(y,z)
%                                         4:  3D line defined on (x,y,z)
%                                                v=f(x,y,z)
%                                                Default is 4
%                       FElines.datapacking --- whether it is block or point
%                                               (deprecated and set to zero
%                                               by default)
%                       FElines.standID ---strandID associated with the
%                                          zone. (Optional).
%                                          If not provided,default is -2
%                                          so that tecplot will assign one
%                                          automatically. Values can be:
%                                          -2 = strand ID to be assigned by tecplot.
%                                          -1 = static strand ID
%                                          0 <= value < 32700 is valid strand ID
%                       FElines.solutiontime --- solution time associated
%                                                with the FEline zone (Optional)
%                                                Default is 0
%                       FElines.zonename --- zone name of each line, if not
%                                       given, then it will be named as 
%                                       FELine1,FELine2,...
%                       FElines.auxname --- auxiliary data name for each 
%                                           line,(cell array of strings
%                                           or single string) (Optional)
%                       FElines.auxval  --- auxiliary data value for each 
%                                           auxname (cell array of strings
%                                           or single string) (Optional)
%                       
%     tdata.FEsurfaces ---an array of finite element 2D surface data of 
%                         triangle or quadrilateral or polygon elements defined
%                         on 2D XY coordinates or 3D XYZ coordinates. Each
%                         surface is treated as zone. FEsurfaces contain
%                         the following information to define the zones:
%                         
%                         FEsurfaces.e2n     ---element connectivity list. For
%                                               each element, it consists of 
%                                               3 (triangular) or 4 (quadrilateral) 
%                                               node numbes, or (polygon) with 
%                                               arbitrary number of points
%                                               (>=3). In case of a
%                                               polygon, e2n is structure
%                                               array rather than an
%                                               integer array. In the
%                                               structure e2n(icell).nodes
%                                               give the 1D array of node
%                                               numbers for each cell
%                                               number icell, where
%                                               1<=icell<=NumElements. 
%                                               For triangular or
%                                               quadrilateral elements, e2n
%                                               has size [NE,3] or [NE,4]
%                                               where NE is number of
%                                               elements. Each row
%                                               (row=1:NE) gives 3 or 4
%                                               node numbers of the
%                                               element.
%
%                         FEsurfaces.f2n     ---face to node connectivity list
%                                               This is not required for FETRIANGLE
%                                               and FEQUADRILATERAL, but can
%                                               be provided to speed up tecplot.
%                                               For FEPOLYGON, this IS required
%                                               as polygons have user-defined 
%                                               number of nodes and edges(face) for 
%                                               each polygonal element.
%                                               f2n is an array of size [NFx2] for
%                                               FETRIANGLE and FEQUADRILATERAL,
%                                               where NF=NE*3 for FETRIANGLE and 
%                                               NF=NE*4 for FEQUADRILATERAL, as
%                                               each triangular element has 3 faces
%                                               (edges) and each quadrilateral elem-
%                                               ent has 4 faces (edges). 
%                                               For FEPOLYGON, f2n is a 1D structure
%                                               array of size [NFx1], and 
%                                               f2n(iface).nodes(1:2) gives the node
%                                               numbers of face number iface.
%                                               1<=iface<=NF. NF is total number of
%                                               faces. 
%                                               ***f2n is not yet implemented***
%
%                         FEsurfaces.varloc  ---location of variable (0-nodal)
%                                               1-cell center (center of element))
%                         FEsurfaces.v       ---value of each variable(nodal or cell)
%                                               2D [nv,Nne]. where nv is number
%                                               of variables in v, Nne=Nn
%                                               if nodal. Nne=Ne if cell-center
%                                               Nn is number of nodes
%                                               Ne is number of elements
% 
%                         FEsurfaces.x      ---x coordinate (1D array)
%                         FEsurfaces.y      ---y coordinate (1D array)
%                         FEsurfaces.z      ---z coordinate (1D array)
%
%                         FEsurfaces.order  ---order of each surface, 
%                                         3 : 2D surface defined on (x,y)
%                                                z=f(x,y),v=f(x,y),z is set
%                                                to zero if not given
%                                         2 : 2D surface defined on (x,z)
%                                                y=f(x,z),v=f(x,z),y is set
%                                                to zero if not given
%                                         1:  2D surface defined on (y,z)
%                                                x=f(y,z),v=f(y,z), x is
%                                                set to zero if not given
%                                         4:  3D surface defined on (x,y,z)
%                                                v=f(x,y,z)
%                                                Default is 4.
%
%                         FEsurfaces.auxname ---auxiliary variable name
%                                        (cell array of strings or single 
%                                        string)
%                         FEsurfaces.auxval  ---auxiliary variable value
%                                        (cell array of strings or single 
%                                        string)
%                         FEsurfaces.datapacking --- is it block or point
%                                                  (deprecated and set to
%                                                  zero by default)
%                         FEsurfaces.strandID ---strandID of the zone. 
%                                                Optional. 
%                                          If not provided,default is -2
%                                          so that tecplot will assign one
%                                          automatically. Values can be:
%                                          -2 = strand ID to be assigned by tecplot.
%                                          -1 = static strand ID
%                                          0 <= value < 32700 is valid strand ID
%                         FEsurfaces.solutiontime --- solution time
%                                                associated with the zone
%                                                Optional, default =0
%                         
%     tdata.FEvolumes  ---an array of finite element 3D volume data with 
%                         tetrahedron elements (4 points per element) or
%                         brick elements (8 points per element) defined on 3D
%                         XYZ coordinates. Each finite element volume is
%                         treated as a zone. FEvolumes contain the
%                         following information to define the zones:
%
%                     
%                         FEvolumes.e2n     ---element connectivity list. For
%                                              each element, it consists of 
%                                              4 (tetrahedron)
%                                              or 8 (brick) node numbers.
%                                              For FEPOLYHEDRON elements,
%                                              e2n is a structure vector,
%                                              with e2n(icell).nodes
%                                              storing node numbers of
%                                              element indexed by icell. 
%
%                         FEvolumes.f2n     ---face to node connectivity list
%                                               This is not required for FETETRAHEDRON
%                                               FEBRICK, but can  be provided to 
%                                               speed up tecplot.
%                                               For FEPOLYHEDRON, this IS required
%                                               as polyhedrons have user-defined 
%                                               number of nodes and faces for 
%                                               each polyhedron element.
%                                               f2n is an array of size [NFx3] for
%                                               FETETRAHEDRON and [NFx4] for 
%                                               FEBRICK, 
%                                               where NF=NE*4 for FETETRAHEDRON and 
%                                               NF=NE*6 for FEBRICK, as
%                                               each tetrahedron element has 4 faces
%                                               and each quadrilateral (brick) elem-
%                                               ent has 6 faces. Each face has 
%                                               3 nodes for FETETRAHEDRON and 4 nodes
%                                               for FEBRICK. 
%                                               For FEPOLYHEDRON, f2n is a 1D structure
%                                               array of size [NFx1], and 
%                                               f2n(iface).nodes(1:NNOF) gives the node
%                                               numbers of face number iface.
%                                               1<=iface<=NF. NF is total number of
%                                               faces. Where NNOF is number of nodes
%                                               on face number iface. 
%                                               ****f2n is not yet implemented*****
%
%                         FEvolumes.varloc  ---location of variable (0-nodal,
%                                              1-cell center)
%                         FEvolumes.v       ---value of each variable(nodal or cell)
%                                              2D [nv,Nne]. where nv is number
%                                              of variables in v, Nne=Nn
%                                              if nodal. Nne=Ne if cell-center
%                                              Nn is number of nodes
%                                              Ne is number of elements
%
%                         FEvolumes.x       ---x coordinate of each node
%                         FEvolumes.y       ---y coordinate of each node
%                         FEvolumes.z       ---z coordinate of each node
%                         FEvolumes.auxname ---auxiliary variable name
%                                        (cell array of strings or single 
%                                        string)
%                         FEvolumes.auxval  ---auxiliary variable value
%                                        (cell array of strings or single 
%                                        string)
%                         FEvolumes.datapacking --- is it block or point
%                                                   (deprecated and set to
%                                                   zero by default)
%                         FEvolumes.strandID --- standID associated with
%                                                the zone (Optional).
%                                          If not provided,default is -2
%                                          so that tecplot will assign one
%                                          automatically. Values can be:
%                                          -2 = strand ID to be assigned by tecplot.
%                                          -1 = static strand ID
%                                          0 <= value < 32700 is valid strand ID
%                         FEvolumes.solutiontime --- solutiontime associated with 
%                                                    the zone. (Optional)
%
%
%      tdata.texts      --- a structure arry of texts that can coexsist with the 
%                           zones. Each element of texts is treated
%                           as a text record (similar to zone). Each text
%                           record can be associated with a zone. texts
%                           contain the following information to define the
%                           text records:
%
%                           texts.str      --- string value of the text
%                                              Must have it.
%                                              Multiple line text
%                                              is given by 'line1 \n
%                                              line2' with "\n" separating
%                                              different lines
%
%                           texts.cs      ---coordinate system of the text
%                                                0= Grid,
%                                                1= Frame,
%                                                2= FrameOffset(not used)
%                                                3= OldWindow(not used)
%                                                4= Grid3D(New to V10)
%                                                Default is 0
%                           texts.scope   ---scope of the text
%                                                0-global, 1-local
%                                                Default is 0.
%
%                           texts.x       -- x coordinate of 
%                                            starting location
%                           texts.y       -- y coordinate of 
%                                            starting location
%                           texts.z       -- z coordinate of starting 
%                                            location (Default to zero)
%                                            Only used for Grid3D coordiante
%                                            system
%                           texts.theta   -- polar angle of starting location
%                                            when x is not given
%                           texts.r       -- polar radius of starting location
%                                            when y is not given
%                           texts.fonttype-- Font of the text
%                                             0-Font_Helvetica,
%                                             1-Font_HelveticaBold,
%                                             2-Font_Greek, 
%                                             3-Font_Math, 
%                                             4-Font_UserDefined, 
%                                             5-Font_Times, 
%                                             6-Font_TimesItalic,
%                                             7-Font_TimesItalicBold, 
%                                             8-Font_TimesBold, 
%                                             9-Font_Courier, 
%                                            10-Font_CourierBold
%                                            Default is 5
%                           texts.charheightunit
%                                         -- unit of font size
%                                               0-grid, 1-frame, 2-Point
%                                               Default is Point
%                           texts.charheight -- height of the characters
%                                               Default is 12. E.g. 12 Point
%                                               if charheightunit is 2
%                           texts.boxtype -- to draw box around text
%                                             0-NOBOX, 1-FILLED, 2-HOLLOW
%                                             Default is 0
%
%                           texts.boxmargin       ---margin of text box
%                                                    Default is 1
%                           texts.boxlinewidth    ---line thickness
%                                                    Default is 1
%                           texts.boxlinecolor    ---outline color of box
%                                                    Default is 0 (black)
%                           texts.boxfillcolor    ---fill color of box
%                                                    Default is 7 (white)
%                           texts.angle   ---angle of texts in degrees
%                                            Default is 0
%               
%                           texts.linespace --- spacing for multi-line
%                                                texts,1-single, 2-double,
%                                                1.5- one and half, etc
%                                                Multiple line text
%                                                is given by "line1 \\n
%                                                line2" with "\\n" separating
%                                                different lines
%                                                Default is 1
%                            
%                           texts.anchor  ---position of the anchor point
%                                            There are nine possible
%                                            anchor positions:
%                                              0=left,       1=center,
%                                              2=right,      3=midleft,
%                                              4=midcenter,  5=midright,
%                                              6=headleft,   7=headcenter,
%                                              8=headright
%                                            Default is 0
%
%                           texts.zone      ---zone that the text will
%                                              be attached to.
%                                              Defalut is 0 (all)
%
%                           texts.color    --- color of the text
%                                              0-Black, 1-R, 2-G,3-B
%                                              4-Cyan 5-Yellow
%                                              6-Purple,7-White
%                                              Defalut is 0
%                           texts.clipping --- clipping method
%                                              0-ClipToAxes
%                                              1-ClipToViewport
%                                              2-ClipToFrame
%                                              Default is 0
%
%       tdata.geometries  --- an array of geometries that can coexist with 
%                           data zones. Each element of geometries array
%                           is treated as a geometry record. Each geometry
%                           record can be associated with a zone.
%                           geometries contain the following information to
%                           define the geometry records:
%
%                           geometries.geomtype ---type of geometry
%                                                  which can be SQUARE,
%                                                  RECTANGLE,CIRCLE,
%                                                  ELLIPSE,LINE,LINE3D
%                                                  0=Line,   1=Rectangle 
%                                                  2=Square, 3=Circle,
%                                                  4=ellipse,5=Line3D
%                                                  Must have. Where
%                                                  Line3D is only
%                                                  available when 
%                                                  geometries.cs is Grid3D
%
%                           geometries.cs  --- coordinate system of 
%                                              geometry, which can be
%                                                0= Grid,
%                                                1= Frame,
%                                                2= FrameOffset(not used)
%                                                3= OldWindow(not used)
%                                                4= Grid3D(New to V10)
%                                                Default is 0
%
%                           geometries.scope --- Scope of the geometry
%                                                0-global, 1-local
%                                                Default is 0.
%
%                           geometries.draworder --- order of drawing
%                                                0-after, 1-before
%                                                Default is 0.
%
%                           geometries.zone  --- zone number to attach
%                                         geometry to. (0 for all zones)
%                                         Default is 0
%
%                           geometries.color --- color of the geometry
%                                                0-Black, 1-Red,  2-Green
%                                                3-Blue,  4-Cyan, 5-Yellow
%                                                6-Purple,7-White
%                                                Default is 0.
%
%                           geometries.fillcolor
%                                            ---fill color of the geometry
%                                               0-Black,1-Red,  2-Green,
%                                               3-Blue, 4-Cyan, 5-Yellow,
%                                               6-Purple,7-White. Default
%                                               is 0.
%
%                           geometries.isfilled
%                                            ---0 = not filled
%                                               1 = filled
%                                               Default is 0
%
%                           geometries.x   --- x coordinate of geometry's
%                                                anchor point 
%                           geometries.y   --- y coordinate of geometry's
%                                                anchor point
%                           geometrrie.z   --- z coordiante of LINE3D's
%                                                anchor point. Default is 0
%                           geometries.theta --- angle of polar line 
%                                                when x  not given
%                                                (used for polar
%                                                coordinate). When both x
%                                                and theta available, x is
%                                                used.
%                           geometries.r    --- radius of polar line
%                                                when y not given
%                                                (used for polar
%                                                coordinate). When both y
%                                                and r available, y is
%                                                used.
%                           geometries.linepattern
%                                           ---   pattern of lines used
%                                                 to draw the geometry,
%                                                 which can be: 
%                                                 0-SOLID,   1-DASHED,
%                                                 2-DOTTED,  3-DASHDOT,
%                                                 4-LONGDASH,
%                                                 6-DASHDOTDOT
%                                                 Default is 0
%
%                           geometries.linethickness --- thickness of 
%                                        drawing line in frame units.
%                                        Default is 0.001, i.e. 0.1%
%
%                           geometries.patternlength --- patternleinght
%                                                  in frame units. 
%                                                  Default is 0.01. 
%
%                           geometries.data  ---data array that describes
%                                               the geometry
%                                         for SQUARE - one number of
%                                                side length
%                                         for RECTANGLE - two numbers of       
%                                                width and height
%                                         for CIRCLE - one number of radius
%
%                                         for ELLIPSE - two numbers of hor-
%                                             izontal axis and vertical 
%                                             axis length
%
%                                         for LINE or LINE3D, the data is a
%                                         colloection of up to 50 polylines
%                                         and each polyline is defined by 
%                                         a number of XY or XYZ coordinates
%                                         The cooridnates are relative
%                                         to the anchor point.
%                                          
%                              %    data(iline).x -- x coordinate vector
%                              %    data(iline).y -- y coordinate vector
%                              %    data(ilien).z -- z coordinate vector
%                                          z is only required for LINE3D
%
%                            geometries.datapacking --point or block
%                                         for POINT, the LINE and LINE3D
%                                         data are written point by point
%                                         for BLOCK, the LINE and LINE3D
%                                         data are written block by block
%                                         (deprecated and set to zero by
%                                         default)
%
%                            geometries.datatype -- datatype (DT) of the 
%                                        point coordinates for LINE and
%                                        LINE3D data. (1-float, 2-double)
%                                        Default 1.
%
%                            geometries.numellipsepts -- number of points
%                                         used to draw circles and ellipses
%                                         for CIRCLE and ELLIPSE types. 
%                                         Default value is 72
%
%                            geometries.clipping -- method of cliping 
%                                                   the geometry.
%                                                   0=ClipToAxes
%                                                   1=ClipToViewport
%                                                   2=ClipToFrame
%                                                   Default is 0
%
%                      Following are only related to LINE geometry
%                            geometries.arrowheadstyle  ---0-plain
%                                                          1-hollow
%                                                          2-filled
%                                                          Default 0
%                            geometries.arrowheadattach ---0-none
%                                                          1-beginning
%                                                          2-end
%                                                          3-both
%                                                          Default 0
%                            geometries.arrowheadsize   ---size of arrow
%                                                          head in frame units
%                                                          Default 0.01
%                            geometries.arrowheadangle  ---angle of arrow head 
%                                                          in degrees. Default 0
%
%
%       tdata.customlabels   --- cell array of custom labels for axis, contour
%                                legend, value or node labels etc.
%
%       tdata.auxdata        --- Meta data structure array that are used to
%                                describe the whole dataset file.
%
%                                auxdata.name  ---name of an aux data
%                                auxdata.value ---value of an aux data
%                                                 (string only).
%
% Input:
%        tdata   --- A structure that contains all kinds of data structure
%                    that tecplot can visualize
%        tecfile --- filename to be created. 
%
% Output: 
%
%        s       ---  0 if successful
%                    -1 if unsuccessful
%
% Example 1:  generate a 3D line in (x,y,z) with temperature 
%             T defined on (x,y,z)
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.lines(1).zonename='myline zone'
%      tdata.lines(1).x=[1,2,3,4,5,6,7,8,9,10];
%      tdata.lines(1).y=[1,2,3,4,5,6,7,8,9,10]+10;
%      tdata.lines(1).z=[1,2,3,4,5,6,7,8,9,10]+1;
%      tdata.lines(1).v(1,:)=[10,20,30,30,20,10,10,10,20,20];
%      mat2tecplot(tdata,'myline3D.plt')
%
% Example 2:  generate a 2D line in (x,y) with temperaature T deined on
%             (x,y). Note there is no z in this case. By default
%             line is defined by points of (x,y) pairs, and 
%             variable T is a function defined on the line.
%      tdata=[];
%      tdata.Nvar=3;
%      tdata.varnames={'x','y','T'};
%      tdata.lines(1).zonename='myline zone';
%      tdata.lines(1).x=[1,2,3,4,5,6,7,8,9,10];
%      tdata.lines(1).y=[1,2,3,4,5,6,7,8,9,10]+10;
%      tdata.lines(1).v(1,:)=[10,20,30,30,20,10,10,10,20,20];
%      mat2tecplot(tdata,'myline2DT_xy.plt')
%
% Example 3:  generate a 2D surface in (x,y) with temperaature T deined on
%             (x,y)
%      tdata=[];
%      tdata.Nvar=4;     %z is set to zero automatically
%                        %even z dos not exit in surfaces(1) below
%                        %it has to be accounted in Nvar and
%                        %varnames
%      tdata.varnames={'x','y','z','T'};
%      tdata.surfaces(1).zonename='mysurface zone';
%      tdata.surfaces(1).x=[1,2,3;1,2,3;1,2,3];    %size 3x3 
%      tdata.surfaces(1).y=[1,1,1;2,2,2;3,3,3];    %size 3x3
%      tdata.surfaces(1).v(1,:,:)=[10,20,30;30,20,10;10,10,20]; 
%      mat2tecplot(tdata,'mysurf2DT_xy.plt')
%
% Example 4:  generate a 2D surface in (x,z) with temperaature T deined on
%             (x,z)
%      tdata=[];
%      tdata.Nvar=4;     %y is set to zero automatically
%                        %even y dos not exit in surfaces(1) below
%                        %it has to be accounted in Nvar and
%                        %varnames
%      tdata.varnames={'x','y','z','T'};
%      tdata.surfaces(1).zonename='mysurface zone';
%      tdata.surfaces(1).x=[1,2,3;1,2,3;1,2,3];    %size 3x3 
%      tdata.surfaces(1).z=[1,1,1;2,2,2;3,3,3];    %size 3x3
%      tdata.surfaces(1).v(1,:,:)=[10,20,30;30,20,10;10,10,20]; 
%      tdata.surfaces(1).order=2;   %order of surface is 2 
%                                   %with (x,z) being coordinate
%                                   %variables
%      mat2tecplot(tdata,'mysurf2DT_xz.plt')
%
%      %Now re-do it with T defined on cell center only
%      tdata.surfaces(1).v=[];
%      tdata.surfaces(1).v(1,:,:)=[10,20;30 10]; %only 2x2 elements
%                                                %[Imax-1,Kmax-1]
%      tdata.surfaces(1).varloc=1;  %specify variable at cell center
%      mat2tecplot(tdata,'mysurf2DT_xz_cellcenter.plt')
%
%
% Example 5:  Combine Example 3 and 4, 2 surfaces, one in (x,y) plane
%             (ordre 3 by default), the other in (x,z) plane (order =2). 
%             And add a third surface in (y,z) plane (order =1).
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.surfaces(1).zonename='mysurface zone1';
%      tdata.surfaces(1).x=[1,2,3;1,2,3;1,2,3];    %size 3x3 
%      tdata.surfaces(1).y=[1,1,1;2,2,2;3,3,3];    %size 3x3
%      tdata.surfaces(1).v(1,:,:)=[10,20,30;30,20,10;10,10,20];
%      tdata.surfaces(1).strandID=1; %you don't have to have strandID
%      tdata.surfaces(2).zonename='mysurface zone2';
%      tdata.surfaces(2).x=[1,2,3;1,2,3;1,2,3];    %size 3x3 
%      tdata.surfaces(2).z=[1,1,1;2,2,2;3,3,3];    %size 3x3
%      tdata.surfaces(2).v(1,:,:)=[10,20,30;30,20,10;10,10,20]; 
%      tdata.surfaces(2).strandID=2;
%      tdata.surfaces(2).order=2;   %order of surface is 2 
%                                   %with (x,z) being coordinate
%                                   %variables
%      tdata.surfaces(3).zonename='mysurface zone3';
%      tdata.surfaces(3).y=[1,2,3;1,2,3;1,2,3];    %size 3x3 
%      tdata.surfaces(3).z=[1,1,1;2,2,2;3,3,3];    %size 3x3
%      tdata.surfaces(3).v(1,:,:)=[10,20,30;30,20,10;10,10,20]; 
%      tdata.surfaces(3).strandID=3;
%      tdata.surfaces(3).order=1;   %order of surface is 2 
%                                   %with (x,z) being coordinate
%                                   %variables
%      mat2tecplot(tdata,'mysurf2DT_xy_xz_yz.plt')
%
% Example 6:  A 3D Finite Element line defined by (x,y,z) and temperature
%             T on the line
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FElines(1).zonename='mysurface zone';
%      tdata.FElines(1).x=[0,1,1,1];  %line made of 4 points (nodes)
%      tdata.FElines(1).y=[0,0,1,1];
%      tdata.FElines(1).z=[0,0,0,1];
%      tdata.FElines(1).v(1,:)=[10,20,30,10];%temperature on 4 nodes
%      tdata.FElines(1).e2n=[1,2;2,3;3,4];  %element connectivity
%                                           %of 3 elements (line
%                                           %segment), each row 
%                                           %gives the two node
%                                           %points of an element
%      mat2tecplot(tdata,'myFEline3D.plt')
%
%  Example 7: a 2D Finite Element surface defiend in (x,y) plane, with 
%             Temperature data on it. T=T(x,y)
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FEsurfaces(1).zonename='my surface zone';
%      tdata.FEsurfaces(1).x=[0,1,0.5]; %triangle of 3 points
%      tdata.FEsurfaces(1).y=[0,0,1];  
%      tdata.FEsurfaces(1).order=3;  %surface defiend on (x,y) coord
%      tdata.FEsurfaces(1).e2n(1,:)=[1,2,3];  %only one element
%                                             %one row of 3
%                                             %node numbers
%      tdata.FEsurfaces(1).v(1,:)=[10,20,30];    %temperature on 3 nodes
%      mat2tecplot(tdata,'myFEsurface_xy.plt')
%
%  Example 8: a 2D Finite Element (triangle) defined in (y,z) plane with 
%             Temperature data on it T=T(y,z)
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FEsurfaces(1).zonename='my surface zone';
%      tdata.FEsurfaces(1).y=[1,2,2.5,3.5,4]; %totally 5 nodes 
%      tdata.FEsurfaces(1).z=[1,3,1,5,1];  
%      tdata.FEsurfaces(1).order=1;  %surface defiend on (y,z) coord
%      tdata.FEsurfaces(1).e2n=[1,2,3;3,2,4;3,5,4];%3 elements(row)
%                                        %each with 3 node numbers
%                                        %(each row has 3 columns)
%      tdata.FEsurfaces(1).v(1,:)=[10,20,30,10,10]; 
%                                        %temperature on 5 nodes
%                                        %only one row (temperatue)
%                                        %(the row has 5 columns)
%      mat2tecplot(tdata,'myFEsurface_yz.plt')
%
%      %now re-do it with T defined on cell centers only:
%      tdata.FEsurfaces(1).v=[];
%      tdata.FEsurfaces(1).v(1,:)=[10,20,30];  %only 3 elements 
%      tdata.FEsurfaces(1).varloc=1;  %specify variable location at cell
%                                     %center
%      mat2tecplot(tdata,'myFEsurface_yz_cellcenter.plt')
%
%  Example 9: a 3D Finite Element brick with 12 nodes, 2 elements, and 
%             temperature T=T(x,y,z)
%         
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FEvolumes(1).zonename='my FE volume zone';
%      tdata.FEvolumes(1).x=[0,3,6,0,3,6,0,3,6,0,3,6]; %totally 12 nodes 
%      tdata.FEvolumes(1).y=[0,0,0,6,6,6,0,0,0,6,6,6];  
%      tdata.FEvolumes(1).z=[0,1,3,3,4,6,8,9,11,11,12,14];
%      tdata.FEvolumes(1).e2n=[1,2,8,7,4,5,11,10;2,3,9,8,5,6,12,11];
%                                               %2 elements(row)
%                                               %each with 8 node numbers
%                                               %(8 columns per row)
%      tdata.FEvolumes(1).v(1,:)=[10,20,30,10,10,20,10,10,20,30,10,10]; 
%                                               %temperature on 12 nodes
%                                               %only one row (temperatue)
%                                               %each row has 12 columns
%      mat2tecplot(tdata,'myFEvolume_brick.plt')
%
%  Exampple 10: Same as Example 9, but using IJK-ordered (cube) zonetype
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.cubes(1).zonename='my IJK volume zone';
%                            %totally 12 nodes [3x2x2]
%      tdata.cubes(1).x=reshape([0,3,6,0,3,6,0,3,6,0,3,6],[3,2,2]);
%      tdata.cubes(1).y=reshape([0,0,0,6,6,6,0,0,0,6,6,6],[3,2,2]);  
%      tdata.cubes(1).z=reshape([0,1,3,3,4,6,8,9,11,11,12,14],[3,2,2]);
%      tdata.cubes(1).v(1,:,:,:)=reshape([10,20,30,10,10,20,10,10,20,30,10,10],[3,2,2]); 
%      mat2tecplot(tdata,'mycube_IJK_order.plt')
%
%  Exampple 11: Same as Example 10, but temperature T is cell-centered
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.cubes(1).zonename='my IJK volume zone';
%                            %totally 12 nodes [3x2x2]
%      tdata.cubes(1).x=reshape([0,3,6,0,3,6,0,3,6,0,3,6],[3,2,2]);
%      tdata.cubes(1).y=reshape([0,0,0,6,6,6,0,0,0,6,6,6],[3,2,2]);  
%      tdata.cubes(1).z=reshape([0,1,3,3,4,6,8,9,11,11,12,14],[3,2,2]);
%      tdata.cubes(1).v(1,:,:,:)=reshape([10,20],[2,1,1]); 
%      tdata.cubes(1).varloc=1;   %variables are provided 
%                                 %at cell centers (center of elements)
%                                 %i.e. v has size [Imax-1,Jmax-1,Kmax-1]
%                                 %==[2,1,1], where Imax=3,Jmax=Kmax=2
%      mat2tecplot(tdata,'mycube_IJK_order_cellcenter.plt')
%
%  Example 12: a 3D Finite Element Tetrhedron with 4 nodes, 1 element, and
%               T=T(x,y,z)
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FEvolumes(1).zonename='my FE volume zone';
%      tdata.FEvolumes(1).x=[0,1,0.5,0.5]; %tetrahedron of 4 points
%      tdata.FEvolumes(1).y=[0,0,1  ,0.5];  
%      tdata.FEvolumes(1).z=[0,0,0  ,1  ];  
%      tdata.FEvolumes(1).e2n(1,:)=[1,2,3,4]; %1 elements(row)
%                                         %each with 4 node numbers
%                                         %(column)
%      tdata.FEvolumes(1).v(1,:)=[10,20,30,10]; %T on 4 nodes
%                                         %only one row (temperatue)
%                                         %each row has 4 columns
%       mat2tecplot(tdata,'myFEvolume_tetrahedron.plt')
%
%  Example 13: Same as Example 12, a 3D Finite Element Tetrhedron with
%              4 nodes, 1 element, and T=T(x,y,z), but T is cell-centered
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','T'};
%      tdata.FEvolumes(1).zonename='my FE volume zone';
%      tdata.FEvolumes(1).x=[0,1,0.5,0.5]; %tetrahedron of 4 points
%      tdata.FEvolumes(1).y=[0,0,1  ,0.5];  
%      tdata.FEvolumes(1).z=[0,0,0  ,1  ];  
%      tdata.FEvolumes(1).e2n(1,:)=[1,2,3,4]; %1 elements(row)
%                                         %each with 4 node numbers
%                                         %(column)
%      tdata.FEvolumes(1).v(1,:)=[15]; %T on 1 element center
%                                      %Only one variable (row)
%                                      %only one element (column)
%      tdata.FEvolumes(1).varloc=1;  %variables are cell-centered
%      mat2tecplot(tdata,'myFEvolume_tetrahedron_cellcenter.plt')
%
%  Example 14: A tetramesh example.
%
%        d = [-1 1];
%        [x,y,z] = meshgrid(d,d,d);  % A cube
%        x = [x(:);0]; %[x,y,z] are corners of a cube plus the center at 0.
%        y = [y(:);0];
%        z = [z(:);0]; 
%        dt = DelaunayTri(x,y,z);
%        Tes = dt(:,:);  %Tes gives connectivity of the tetrahedrons
%        X = [x(:) y(:) z(:)];
%        tetramesh(Tes,X);camorbit(20,0)
%        tdata=[];
%        tdata.Nvar=4;
%        tdata.varnames={'x','y','z','T'};
%        tdata.FEvolumes(1).zonename='my tetra mesh zone';
%        tdata.FEvolumes(1).x=x';
%        tdata.FEvolumes(1).y=y';
%        tdata.FEvolumes(1).z=z';
%        %temperature
%        tdata.FEvolumes(1).v(1,:)=linspace(10,30,length(z));
%        tdata.FEvolumes(1).e2n=Tes;
%        mat2tecplot(tdata,'myFEvolume_tetramesh.plt')
%    
%  Example 15: show puget sound FVCOM model grid in 2D xy plane with depth h
%              defined on nodes.
%
%      load psm_grid_v2.mat;  % this data is prepared and conatin following
%
%      %    e2n             13941x5    --- element to node connectivity (:,2:4)    
%      %    h_n              9013x1    --- depth of all nodes
%      %    lld_n            9013x4    --- lat, lon and depth on nodes (:,2:4)
%      %    xy_e            13941x2    --- x, y of element centers (:,1:2)
%      %    xyd_n            9013x4    --- x, y and depth on nodes (:,2:4)
%
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','h'};
%      tdata.FEsurfaces(1).zonename='Puget Sound Model 2.0 grid and depth';
%      tdata.FEsurfaces(1).x=xyd_n(:,2);  %2nd column in xyd_n is x   
%      tdata.FEsurfaces(1).y=xyd_n(:,3);  %3rd                    y
%      tdata.FEsurfaces(1).order=3;       %surface defiend on (x,y) coord
%      tdata.FEsurfaces(1).e2n=e2n(:,2:4);%each row(element) has 3 node numbers
%      tdata.FEsurfaces(1).v(1,:)=h_n';   %depth h of all nodes
%      mat2tecplot(tdata,'psm_grid_v2.plt')
%  
%  Example 16: Re-do Example 15, however, show the grid as 3D with triangular cylinders
%              vertically. Tecplot does not support triangular cylinder elements (prism) directly,
%              but we can fool it around using bricks with two nodes on top side of a
%              brick collapsing into one node, and do the same for bottom side.
%
%      load psm_grid_v2.mat;  % this data is pre-prepared and contains the following 2D grid
%
%      %    e2n             13941x5    --- element to node connectivity (:,2:4)
%      %    h_n              9013x1    --- depth of all nodes
%      %    lld_n            9013x4    --- lat, lon and depth on nodes (:,2:4)
%      %    xy_e            13941x2    --- x, y of element centers (:,1:2)
%      %    xyd_n            9013x4    --- x, y and depth on nodes (:,2:4)
%
%      %make up 3D element to node connectivities
%      NumNodes2D=length(h_n);   %number of nodes in 2D grid
%      NumNodes3D=NumNodes2D*2;  %two layers as a simple example
%      e2n2d=[e2n(:,2:4) e2n(:,4)];   %repeat last 2D node of each element 
%                                     %(4 nodes to make top surface of brick)
%      e2n3d=[e2n2d e2n2d+NumNodes2D];%first 4 are node numbers on surface of brick
%                                     %next 4 are node numbers on bottom of brick
%                                     %(which offsets 2D node numbers by NumNodes2D) 
%      h3d=[0*h_n;h_n];  %first depth on all surface nodes (zero)
%                        %followed by depth on all bottom nodes
%      x2d=xyd_n(:,2);  %2D x coordinate of all 2D nodes
%      y2d=xyd_n(:,3);  %2D y coordinate of all 2D nodes
%      x3d=[x2d;x2d];   %repeat to get two vertical layers (x of 3D nodes)
%      y3d=[y2d;y2d];
%      z3d=-h3d*100000;   %z coodinate = -h, also convert to 10 micrometer units 
%                         %otherwise, tecplot aspect ratio will be strange
%      tdata=[];
%      tdata.Nvar=4;
%      tdata.varnames={'x','y','z','h'};
%      tdata.FEvolumes(1).zonename='Puget Sound Model 2.0 3D grid and height';
%      tdata.FEvolumes(1).x=x3d; 
%      tdata.FEvolumes(1).y=y3d;
%      tdata.FEvolumes(1).z=z3d;
%      tdata.FEvolumes(1).e2n=e2n3d;
%      tdata.FEvolumes(1).v(1,:)=h3d'; 
%      mat2tecplot(tdata,'psm_grid_v2_3D_2layers.plt')
%
%  Example 17: A few geometries and texts in frame coordinates
%
%      tdata=[];
%      tdata.Nvar=2;
%      tdata.geometries(1).geomtype=0; %polylines
%      tdata.geometries(1).cs=1;       %in frame coordinates
%      tdata.geometries(1).x=0.5;      %anchoring in the middle of the frame
%      tdata.geometries(1).y=0.5;      
%      tdata.geometries(1).color=1;    %red color
%      tdata.geometries(1).data(1).x=[0,0.1,0.1,0];  %line coordinates relative to anchor
%      tdata.geometries(1).data(1).y=[0.4,0.3,0.2,0.1];
%      tdata.texts(1).cs=1;
%      tdata.texts(1).x=0.5;
%      tdata.texts(1).y=0.8;
%      tdata.texts(1).str='This is a polyline';
%
%      tdata.geometries(2).geomtype=1;  %rectangle
%      tdata.geometries(2).cs=1;       %in frame coordinates
%      tdata.geometries(2).x=0.2;      %anchoring in lower left corner
%      tdata.geometries(2).y=0.2;
%      tdata.geometries(2).color=0;    %black color
%      tdata.geometries(2).data=[0.2,0.3];  %width and heigh of the rectangle
%      tdata.texts(2).cs=1;
%      tdata.texts(2).x=0.2;
%      tdata.texts(2).y=0.2;
%      tdata.texts(2).str='This is a \n rectangle';  %text in two lines
%
%      tdata.geometries(3).geomtype=3; %circle
%      tdata.geometries(3).cs=1;       %in frame coordinates
%      tdata.geometries(3).x=0.8;      %anchoring in top right corner
%      tdata.geometries(3).y=0.8;
%      tdata.geometries(3).color=5;    %yellow
%      tdata.geometries(3).linethickness=0.01; %1 percent frame unit
%      tdata.geometries(3).data=[0.3]; %radius of the circle
%  
%      tdata.geometries(4).geomtype=2; %square
%      tdata.geometries(4).cs=1;       %in frame coordinates
%      tdata.geometries(4).x=0.6;      %anchoring lower right
%      tdata.geometries(4).y=0.2;
%      tdata.geometries(4).color=4;    %cyan line
%      tdata.geometries(4).linethickness=0.01; %1 percent frame unit
%      tdata.geometries(4).data=[0.3]; %width of the square
%      tdata.geometries(4).isfilled=1; %fill the square
%      tdata.geometries(4).fillcolor=6 %fill with purple 6
%
%      tdata.geometries(5).geomtype=4; %ellipse
%      tdata.geometries(5).cs=1;       %in frame coordinates
%      tdata.geometries(5).x=0.6;      %anchoring lower right
%      tdata.geometries(5).y=0.2;
%      tdata.geometries(5).color=7;     %white line
%      tdata.geometries(5).linethickness=0.01; %1 percent frame unit
%      tdata.geometries(5).data=[0.2,0.1]; %x and y axis length of ellipse
%      tdata.geometries(5).isfilled=1;  %fill the square
%      tdata.geometries(5).fillcolor=3; %fill with blue
%
%      mat2tecplot(tdata,'mygeometry.plt');   
% 
% Wen Long, Seattle, 06/15/2012
%
%

if nargin <2   %if number of arguments less than 2
               %simply quit and give help
    help mat2tecplot;
    s=-1;
    return;
end

%
%PrecisoinByCountOfInt32=2;   % float64 (2 -double precision
%                             %          1 -float
%                             %          2 -double
%                             %          3 -longInt
%                             %          4 -shortInt
%                             %          5 -Byte
%                             %          6 -Bit )
%

output_file_name=tecfile; 

fid_out = fopen(output_file_name,'w');  %should  be 'w' 'b' for big endian
if fid_out <0
    disp('Fail to open the file to write, check permission and path format \n');
end

%I . Header Section

%  Header section i. magic number 
%      Header section starts with Magic number (i.e. Version number)
%
%         +-----------+
%         | "#!TDV112"|       8 Bytes, exact characters "#!TDV112".
%         +-----------+       Version number follows the "V" and
%                             consumes the next 3 characters (for
%                             example: "V75 ", "V101").


magic_number = '#!TDV112';
char_hold = char(magic_number);
l_max = max(size(char_hold));
for ii =1:1:l_max
    fwrite(fid_out,char_hold(ii),'char');
end

% Header section ii. Integer value of 1.
%
%  +-----------+
%  | INT32     |       This is used to determine the byte order
%  +-----------+       of the reader relative to the writer.
%
%
dummy_int32 = 1; %one for little endian?
fwrite(fid_out,dummy_int32,'int32');

% Header section iii. Title and variable names.
%
% 	      +-----------+ 
% 	      | INT32     |      FileType: 0 = FULL, 
% 	      +-----------+                1 = GRID, 
% 				                       2 = SOLUTION 
%         +-----------+
%         | INT32*N   |       The TITLE. (See note 1.)
%         +-----------+
%         +-----------+
%         | INT32     |       Number of variables (NumVar) in the datafile.
%         +-----------+
%         +-----------+
%         | INT32*N   |       Variable names.  N =  L[1] + L[2] + .... L[NumVar]
%         +-----------+       where:
%                                    L[i] = length of the ith variable name + 1
%                                           (for the terminating 0 value).
%                                           (See note 1.)

% filetype
dummy_int32 = 0;                      %
fwrite(fid_out,dummy_int32,'int32');

tdatanames=fieldnames(tdata);  %find all the field names under structure 
                               %tdata

have_title=~isempty(find(strcmp(tdatanames,'title')==1));  %true if have title

if(~have_title)                  %set default title to 'Tecplot Data File' if
   tdata.title='tecplot data' ;  % tdata.title does not exist
end

   %  write title and terminate with a null char
plt_write_string(fid_out, tdata.title);

have_Nvar=~isempty(find(strcmp(tdatanames,'Nvar')==1));  %check if have 'Nvar'

   % Number of variables
if(~have_Nvar)  %check if number of variables given
                              %if not given then calculate
                              %from tdata.varnames
    have_varnames=~isempty(find(strcmp(tdatanames,'varnames')==1));
    if(have_varnames)
        tdata.Nvar=length(tdata.varnames); %
    else %if tdata.varnames is not there, then give error and quit
        s=-1;
        display('Error: must give tdata.Nvar or tdata.varnames or both'); 
        return;
    end
    
else
   have_varnames=~isempty(find(strcmp(tdatanames,'varnames')==1));
   if(have_varnames) %
      if(length(tdata.varnames)<tdata.Nvar) %set the remaining variable
                                            %names, if names not enough
         for ivar=length(tdata.varnames):tdata.Nvar
             tdata.varnames{ivar}=['var' int2str(ivar)];
         end
      end
   else
      %set all variable names to var1,var2,...,varNvar
      for ivar=1:tdata.Nvar
          tdata.varnames{ivar}=['var' int2str(ivar)];
      end
   end
end
    
var_count=tdata.Nvar;    
dummy_int32 = var_count;
fwrite(fid_out,dummy_int32,'int32');

%variable names 
%plt_write_string(fid_out, x_name);
%plt_write_string(fid_out, y_name);
for ivar=1:tdata.Nvar 
      plt_write_string(fid_out, tdata.varnames{ivar});
end 
   
   % additional var name list
   %for i=1:length(matvarlist)
   %    plt_write_string(fid_out, matvarlist{i});
   %end

%
% Header section iv, definition of Zones 
%
%   Here is how zones are defined:
%
%         +-----------+
%         | FLOAT32   |       Zone marker. Value = 299.0
%         +-----------+
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).
%         +-----------+
%         | INT32     |       ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+                3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                                      6=POLYGON,7=POLYHEDRON
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+
%
%         +-----------+
%         | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
%         +-----------+       located at the nodes.  1 = Specify
%
%         if "specify var location" == 1
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%
%
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.
%
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
%
%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+
%
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
% 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
% 		        +-----------+ 
% 		        | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32     |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+
%
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      
%
%
%Loop through all liness and treat each line as a zone and define zones
%

have_lines=~isempty(find(strcmp(tdatanames,'lines')==1));  %check if have lines
                
%get number of lines in tdata
if(have_lines)  %make sure have lines
    if(isstruct(tdata.lines))  %make sure tdata.lines is structure array
                               %or at least structure 
       Nlines =length(tdata.lines);
    else
       Nlines = 0;
    end
else
    Nlines=0;   
end

for iline=1:Nlines
  
linefields=fieldnames(tdata.lines(iline));  %find fields in
                                            %lines(iline) structure
% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name 'data zone' and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

have_zonename=~isempty(find(strcmp(linefields,'zonename')==1));
if(have_zonename && isempty(tdata.lines(iline).zonename))
   have_zonename=0;
end
if(~have_zonename)  %if no zone name, give one by default
                    %as 'Line1','Line2' etc, depending on index iline
   zone_name = ['Line' int2str(iline)];  
else
   if(   ischar(tdata.lines(iline).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.lines(iline).zonename))) %and not an empty string
       zone_name=tdata.lines(iline).zonename; 
   else
       %give warning and make up zone name
       warning( ['lines number ' int2str(iline) 'zonename is NOT a string!']);
       zone_name=['Line' int2str(iline)];  %make up names by default
       warning( ['set to default value as ' zone_name]);
   end
end

%output zone_name 
plt_write_string(fid_out, zone_name);     %Here N = length(zone_name)+1

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

dummy_int32 = -1;  %always give -1, assuming no parent zones
fwrite(fid_out,dummy_int32,'int32');

% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

%check if this line has strandID
have_strandID=~isempty(find(strcmp(linefields,'strandID')==1)); 
if(have_strandID && isempty(tdata.lines(iline).strandID))
    have_strandID=0;
end
if(~have_strandID) %if not exist, then give zero by default
   strand_ID= -2;
else
    if(isinteger(tdata.lines(iline).strandID))  %make sure it is an integer
        strand_ID=tdata.lines(iline).strandID;
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.lines(iline).strandID);
        catch 
            strand_ID=-2; %default to -2
        end
        warning(['strandID of line ' int2str(iline) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end
fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

%check if this line has solutiontime
have_solutiontime=~isempty(find(strcmp(linefields,'solutiontime')==1));
if(have_solutiontime && isempty(tdata.lines(iline).solutiontime))
    have_solutiontime=0;
end
if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(  isnumeric(tdata.lines(iline).solutiontime) ...
       &&  isfinite(tdata.lines(iline).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.lines(iline).solutiontime;
    else
        warning(['solutiontime of line ' int2str(iline) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');

%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

%this is default to be deterimned by tecplot
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE,
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON,
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON
%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%
%

%for lines, zone_type is 0, i.e. ordered
zone_type = 0;  %ordered 
fwrite(fid_out,zone_type,'int32');

%
%Data packing
% 0- Block
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

% % check if there is datapacking 
% have_datapacking=~isempty(find(strcmp(linefields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.lines(iline).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.lines(iline).datapacking) ...
%       &&  isfinite(tdata.lines(iline).datapacking) ...
%       &&  (  tdata.lines(iline).datapacking ==0 ...
%           || tdata.lines(iline).datapacking ==1))
%         data_packing=tdata.lines(iline).datapacking;
%    else
%        warning(['datapacking of line ' int2str(iline) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end

data_packing=0;

% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------

%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  not specifying
%   1 ----  Specify
%  

%for lines, all data must be on nodes.

 var_specified=0;  %do not specify
 fwrite(fid_out,var_specified,'int32');

 if(var_specified==1)
%
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)
%    %
%    %give location for each variable in this zone
%    %
%    %for iv=1:NV %NV is number of variables
%        %if variable iv is at node
%           var_location = 0; %Node
%           %fwrite(fid_out,0,'int32'); %if
%        %elseif variable iv is cell centered
%           var_location = 0; %cell centered
%           %fwrite(fid_out,1,'int32');
%        % %end
%    %end
 end
 

%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%

% Here we set it to false for ordered data 
% false = 0 (int32), true = 1 (int32)
%

%for lines, this is set to zero
dummy_int32 = 0; %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,dummy_int32,'int32');

% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

%for lines, this is set to zero (no faces)

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

% if (NoOfUserDefinedNeighbourConn ~=0)  
%     
% %         if "number of miscellaneous user defined face neighbor connections" != 0
% %           +-----------+
% %           | INT32     |     User defined face neighbor mode
% %           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
% %                             3=Global 1-to-many)
% 
%      %Pls modify based on data
%      user_defined_face_neighbor_mode= 0; %Local 1-to-1
%      fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
%      
%      if(zone_type>0)  %for non-oredered zones only
% %           if FE Zone:
% %             +-----------+
% %             | INT32     |     Indicates if the finite element face neighbors are
% %             +-----------+     completely specified by the miscellaneous face neighbors
% %                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
% %                               auto assignment of face neighbors otherwise all faces not
% %                               specified are considered boundaries. If no, then Tecplot will
% %                               perform auto-assignment of the face neighbors unless the
% %                               raw face neighbor array was supplied. This option is not
% %                               valid for ORDERED zones.
%         
%         %pls modify this value based on your data
%         fe_face_complete_by_misc_face_neighbors=1;
%         
%         fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
%      end
%     
%     
% end

if(zone_type==0)  %for ordered zone

%give maximum indices (size) of the dimensions
    
%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   %for line in 2D, it is (x,y),(x,z), or (y,z) 
   %x-y-z, i-j-k, for line in (x,y) plane, z value is given as zero
   %              for line in (x,z) plane, y value is given as zero
   %              for line in (y,z) plane, x value is given as zero
   %
   %for line in 3D, it is (x,y,z), x, y, z values are given
   %
   
   % 
   % ordered zone, 3 int32s for num points i -j -k
   %
   
   %check if it is (x,y) or (x,z) or (y,z)
   %if x,y,z all exist, it is a 3d line, give error and quit
   %
   have_x=~isempty(find(strcmp(linefields,'x')==1));
   have_y=~isempty(find(strcmp(linefields,'y')==1));
   have_z=~isempty(find(strcmp(linefields,'z')==1));

   if(have_x && isempty(tdata.lines(iline).x))
      have_x=0;
   end
   if(have_y && isempty(tdata.lines(iline).y))
      have_y=0;
   end
   if(have_z && isempty(tdata.lines(iline).z))
      have_z=0;
   end

   if(have_x && have_y && have_z)
      Ixmax=length(tdata.lines(iline).x);       
      Iymax=length(tdata.lines(iline).y);
      Izmax=length(tdata.lines(iline).z);
      Imax=min([Ixmax, Iymax,Izmax]);
      if(Ixmax ~= Imax || Iymax~=Imax || Izmax~=Imax)
             warning(['Warning: line ' int2str(iline) 'x, y and z length not equal!']);
             warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
      end
      x_data=tdata.lines(iline).x(1:Imax);
      y_data=tdata.lines(iline).y(1:Imax);
      z_data=tdata.lines(iline).z(1:Imax);
      
   else
   
	   if(~have_z)
	      %make sure size of x equal to size of y for (x,y) line (z not given)
	      %set default z value to zero
	      if(~have_x || ~have_y)
		 s=-1; 
		 display(['Error: line ' int2str(iline) ' x or y value must be given if you do not give z value']);
		 return;
	      else
		 Ixmax=length(tdata.lines(iline).x);
		 Iymax=length(tdata.lines(iline).y);
		 Imax=min(Ixmax,Iymax); %take minimum value
		 if(Ixmax ~= Iymax)
		     warning(['Warning: line ' int2str(iline) ' x and y length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 x_data=tdata.lines(iline).x(1:Imax);
		 y_data=tdata.lines(iline).y(1:Imax);
		 z_data=zeros(size(x_data)) ; %default to zero for z when it does not exist
	      end      
	   end
	   if(~have_y)
	      %make sure size of x equal to size of z for (x,z) line (y not given)
	      %set default y value to zero
	      if(~have_x || ~have_z)
		 s=-1; 
		 display(['Error: line ' int2str(iline) ' x or z value must be given if you do not give y value']);
		 return;
	      else
		 Ixmax=length(tdata.lines(iline).x);
		 Izmax=length(tdata.lines(iline).z);
		 Imax=min(Ixmax,Izmax); %take minimum value
		 if(Ixmax ~= Izmax)
		     warning(['Warning: line ' int2str(iline) ' x and z length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 x_data=tdata.lines(iline).x(1:Imax);
		 z_data=tdata.lines(iline).z(1:Imax);
		 y_data=zeros(size(x_data)) ; %default to zero for y when it does not exist
	      end
	   end
	   if(~have_x)
	      %make sure size of y equal to size of z for (y,z) line (x not given)
	      %set default x value to zero
	      if(~have_y || ~have_z)
		 s=-1;           
		 display(['Error: line ' int2str(iline) ' y or z value must be given if you do not give x value']);
		 return;
	      else
		 Izmax=length(tdata.lines(iline).z);
		 Iymax=length(tdata.lines(iline).y);
		 Imax=min(Izmax,Iymax); %take minimum value
		 if(Izmax ~= Iymax)
		     warning(['Warning: line ' int2str(iline) ' y and z length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 z_data=tdata.lines(iline).z(1:Imax);
		 y_data=tdata.lines(iline).y(1:Imax);
		 x_data=zeros(size(z_data)) ; %default to zero for x when it does not exist
	      end
	   end
   end
   
   %
   %first dimension is number of points on the line
   %it is always the length of x or y or z
   %
   
   %Imax  is calculated above based on number of points in x, y, z
   
   Jmax = 1; %J is set to 1
   Kmax = 1; %K is set to 1
   fwrite(fid_out,Imax,'int32');
   fwrite(fid_out,Jmax,'int32');
   fwrite(fid_out,Kmax,'int32');

   %this assumes that you always have "x","y","z" as variables in 
   %the output file, other variables might be vectors or scalars
   %defined on the points of the lines, such as pressure, temperature,
   %salinity, velocity etc. It really depends on what variables you want to
   %have. If they are given, then we output the first Imax values of them.
   %If they are not given, then we output them as zero
   
else  %for non-ordered (FE) zones
   
    %for lines, this section is ignored
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
% 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary connections. 
% 			    +-----------+
%
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       
%      % %uncomment and complete here
% 
%       %NumPts=  ; %give number of nodes
%       %fwrite(fid_out,NumPts,'int32');
         if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      

%          %uncomment the following and complete
%          %
%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');
         end      
%       
%       %NumElements= ; %number of elements
%       %fwrite(fid_out,NumElements,'int32');
%       
%       ICellDim=0;
%       JCellDim=0;
%       KCellDim=0;
%       fwrite(fid_out,ICellDim,'int32');
%       fwrite(fid_out,JCellDim,'int32');
%       fwrite(fid_out,KCellDim,'int32');
    
end


% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%
% Here we only dealw ith texts for 
%       tdata.lines(iline).auxname and tdata.lines(iline).auxval 
%

%check if we have aux data for this line, true if both auxname and 
%auxval exist and are 'cell arrays or string' and neither of them is empty
have_auxdata=(   ~isempty(find(strcmp(linefields,'auxname')==1))  ...
              && ~isempty(find(strcmp(linefields,'auxval')==1)));   
          
if(have_auxdata)

   aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
   aux_names=tdata.lines(iline).auxname;  %cell array or string
   aux_values=tdata.lines(iline).auxval; %cell array or string
   %
   %check if they are cell arrays (1D) or just a single string
   %
   if(iscellstr(aux_names) && iscellstr(aux_values))
      %or check if they are string (not array), just single string value
      N_aux=min(length(aux_names),length(aux_values));
                             %number of aux is obtained from aux_names
                             %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values)> N_aux)
          %give warnging if they do not match in size
          warning(['Warning: line ' int2str(iline) ' auxname and auxval have inconsistent sizes']);
          warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
   else
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
   end
   
   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: line ' int2str(iline) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this line']);
       aux_data=0;
   end
   
else
   aux_data=0;  %no aux data
end


%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.
if(aux_data==0)
    fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
            
          fwrite(fid_out,aux_data,'int32');  %repeat for each aux variable
          
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
              aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value
          
    end
end 

end %end of iline loop


%Loop through all surfaces, and treat each surface as a zone

%----Surface Data---
have_surfaces=~isempty(find(strcmp(tdatanames,'surfaces')==1));  %check if have surfaces
%get number of surfaces in tdata
if(have_surfaces)  %make sure have surfaces
    if(isstruct(tdata.surfaces))  %make sure tdata.surfaces is structure array
                                  %or at least structure 
       Nsurfs =length(tdata.surfaces);
    else
       Nsurfs = 0;
    end
else
    Nsurfs=0;   
end

for isurf=1:Nsurfs


surffields=fieldnames(tdata.surfaces(isurf));  %obtain field names
                                               %in surfaces(isurf)

% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

%check if have zonenames for this surface

have_zonename=~isempty(find(strcmp(surffields,'zonename')==1));

if(have_zonename && isempty(tdata.surfaces(isurf).zonename))
   have_zonename=0;
end

if(~have_zonename)  %if no zone name, give one by default
                    %as 'Surface1','Surface2' etc, depending on index isurf
   zone_name = ['Surface' int2str(isurf)];  
else
   if(   ischar(tdata.surfaces(isurf).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.surfaces(isurf).zonename))) %and not an empty string
       zone_name=tdata.surfaces(isurf).zonename; 
   else
       %give warning and make up zone name
       warning(['surfaces number ' int2str(isurf) 'zonename is NOT a string!']);
       zone_name=['Surface' int2str(isurf)];  %make up names by default
       warning(['set to default value as ' zone_name]);
   end
end

%output zone_name
                          
plt_write_string(fid_out, zone_name);

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

%check if have strandID for this surface zone
have_strandID=~isempty(find(strcmp(surffields,'strandID')==1)); 
if(have_strandID && isempty(tdata.surfaces(isurf).strandID))
    have_strandID=0;
end

if(~have_strandID) %if not exist, then give zero by default
   strand_ID= -2;
else
    if(isinteger(tdata.surfaces(isurf).strandID))  %make sure it is an integer
        strand_ID=tdata.surfaces(isurf).strandID;
        
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.surfaces(isurf).strandID);
        catch 
            strand_ID=-2; %default to  -2
        end
        warning(['strandID of surface ' int2str(isurf) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end
fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

%check if have solutiontime

have_solutiontime=~isempty(find(strcmp(surffields,'solutiontime')==1));
if(have_solutiontime  && isempty(tdata.surfaces(isurf).solutiontime))
    have_solutiontime=0;
end

if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(  isnumeric(tdata.surfaces(isurf).solutiontime) ...
       &&  isfinite(tdata.surfaces(isurf).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.surfaces(isurf).solutiontime;
    else
        warning(['solutiontime of surface ' int2str(isurf) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');


%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

%this default to -1 and to be determined by tecplot
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE,
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON,
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON
%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!

%for surfaces zone_type is 0, i.e. ordered
zone_type = 0;  %ordered
fwrite(fid_out,zone_type,'int32');


%
%Data packing
% 0- Block
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

% % check if there is datapacking 
% have_datapacking=~isempty(find(strcmp(surffields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.surfaces(isurf).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.surfaces(isurf).datapacking) ...
%       &&   isfinite(tdata.surfaces(isurf).datapacking) ...
%       &&  (  tdata.surfaces(isurf).datapacking ==0 ...
%           || tdata.surfaces(isurf).datapacking ==1))
%         data_packing=tdata.surfaces(isurf).datapacking;
%    else
%        warning(['datapacking of surface ' int2str(isurf) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end
data_packing=0;

% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------


%
%find which order the surface data is organized.
%
%If IJ-ordered, then x and y must be
%given nodal, and have size (Imax*Jmax),Kmax must be one. 
%z will be ignored if not supplied. If z is supplied 
%it will be listed as the third variable

%If IK-ordered, then x, z must be given as nodal, and have size
%(IMAX*KMAX), JMAX must be one. y will be ignored if not supplied
%If y is supplied, it will be listed as the second variable
%
%IF JK-ordered, then y,z must be given as nodal, and have size (JMAX*KMAX),
%IMAX must be one. x will be ignored if not supplied. If x is supplied
%it will be listed as the first variable
%
%check if order is available
have_order=~isempty(find(strcmp(surffields,'order')==1));
if(have_order && isempty(tdata.surfaces(isurf).order))
    have_order=0;
end

if(have_order)
    if( isnumeric(tdata.surfaces(isurf).order) ...
       &&(  tdata.surfaces(isurf).order==1  ... 
         || tdata.surfaces(isurf).order==2  ...
         || tdata.surfaces(isurf).order==3) )
         surface_order=tdata.surfaces(isurf).order;
    else
        warning(['Oops, surface ' int2str(isurf) ' order value invalid']);
        warning(['set to default value of 3 (IJ-ordered)']);
        surface_order=3;
    end
else
    %default to have order =3 (IJ-ordered data)
    surface_order=3;
    warning(['Oops, surface ' int2str(isurf) ' order not given']);
    warning(['set to default value of 3 (IJ-ordered),and you must have x, y']);
end
have_x=~isempty(find(strcmp(surffields,'x')==1));
have_y=~isempty(find(strcmp(surffields,'y')==1));
have_z=~isempty(find(strcmp(surffields,'z')==1));
if(have_x && isempty(tdata.surfaces(isurf).x))
   have_x=0;
end
if(have_y && isempty(tdata.surfaces(isurf).y))
   have_y=0;
end
if(have_z && isempty(tdata.surfaces(isurf).z))
   have_z=0;
end
	
%
%Check variable size to make sure coordinate variables are nodal and have same size
%
switch(surface_order)
    case 3
        %IJ-ordered, x, y must have same size [Imax, Jmax]
        if(have_x && have_y)  %must have x and y
            Imax_x=size(tdata.surfaces(isurf).x,1);
            Jmax_x=size(tdata.surfaces(isurf).x,2);
            Imax_y=size(tdata.surfaces(isurf).y,1);
            Jmax_y=size(tdata.surfaces(isurf).y,2);
            if(Imax_x==Imax_y && Jmax_x==Jmax_y)
                Imax=Imax_x;
                Jmax=Jmax_x;
                Kmax=1;
            else %oops size do not match
                s=-1;                
                display(['Error: surface ' int2str(isurf) ' does not have matching x and y dimensions!!']);
                return;
            end
        else
            s=-1;
            error(['Error: surface ' int2str(isurf) ' does not have x or y ']);
            return;
        end
    case 2
        %IK-ordered, x, z must have same size [Imax, Kmax]
        if(have_x && have_z)  %must have x and y
            Imax_x=size(tdata.surfaces(isurf).x,1);
            Kmax_x=size(tdata.surfaces(isurf).x,2);
            Imax_z=size(tdata.surfaces(isurf).z,1);
            Kmax_z=size(tdata.surfaces(isurf).z,2);
            if(Imax_x==Imax_z && Kmax_x==Kmax_z)
                Imax=Imax_x;
                Kmax=Kmax_x;
                Jmax=1;
            else %oops size do not match
                s=-1;
                display(['Error: surface ' int2str(isurf) ' does not have matching x and z dimensions!!']);
                return;
            end
        else
            s=-1;
            display(['Error: surface ' int2str(isurf) ' does not have x or z ']);
            return;
        end
    case 1
        %JK ordered, y, z must have same size (Jmax, Kmax]
        if(have_y && have_z)  %must have y and z
            Jmax_y=size(tdata.surfaces(isurf).y,1);
            Kmax_y=size(tdata.surfaces(isurf).y,2);
            Jmax_z=size(tdata.surfaces(isurf).z,1);
            Kmax_z=size(tdata.surfaces(isurf).z,2);
            if(Jmax_y==Jmax_z && Kmax_y==Kmax_z)
                Jmax=Jmax_y;
                Kmax=Kmax_y;
                Imax=1;
            else %oops size do not match
                s=-1;
                display(['Error: surface ' int2str(isurf) ' does not have matching y and z dimensions!!']);
                return;
            end
        else
            s=-1;
            display(['Error: surface ' int2str(isurf) ' does not have y or z ']);
            return;
        end        
end
%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  Not specifying
%   1 ----  Specify
%  

%
%Wen Long: note that when data_packing=1, i.e. point (IsBlock=0), var
%location must be nodal!! That is Cell Centered variables can only occur
%when data_packing is block. Cell Centered variables CAN'T happen when
%data_packing method is point. For Nodal variables, they can occur for both
%point packing and block packing methods.
%

%check if varloc is available

have_varloc=~isempty(find(strcmp(surffields,'varloc')==1));
if(have_varloc && isempty(tdata.surfaces(isurf).varloc))
   have_varloc=0; 
end
if(~have_varloc) %give default value as 0 (nodal)
   var_loc=0;
   var_specified=0;
else
   if(    isfinite(tdata.surfaces(isurf).varloc) ...
      &&  (  tdata.surfaces(isurf).varloc ==0 ...
          || tdata.surfaces(isurf).varloc ==1))
        var_loc=tdata.surfaces(isurf).varloc;
        var_specified=1;
        
        %Wen Long, make sure when var_specified==1
        %data_packing is 0 (block) 
        %rule out conflict conditions (data_packing=1 (point) and var_loc=
        %1 (cell-center) cannot co-exist)
        %That is to say when var_loc=1, data_packing must be zero (block)
        %
        if(data_packing==1 && var_loc==1)
            s=-1;
            display(['Error: datapacking of surface ' int2str(isurf) ' (1 point) conflicts with varloc']);
            return;
        end
   else
       warning(['var location of surface ' int2str(isurf) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
       var_loc=0;
       var_specified=0;
       warning(['set default value as 0 (nodal)']);
   end
end

 %output var_specified
 fwrite(fid_out,var_specified,'int32');
 
 if(var_specified==1)
     
%
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)

%Current method assumes all variables have the same type of locations
%
% Assume NV=Nn+Nc, where Nn is number of nodal variables
%                        Nc is number of non-nodal variables
% For surface, we must have Nn>=2, 
%
% We treat the surface based on how the surface is defined:
%
% for surface in (x,y) plane, e.g. P(x,y), need x, y be nodal -- IJ-ordered
%                                          z can be given as nodal or not.
%                                          z=f(x,y)
%                                          If z is constant, it is a 2D
%                                          surface. If z is not constant
%                                          z can be used to interpret the
%                                          surface as a 3D surface
%   
% for surface in (x,z) plane, e.g. P(x,z), need x, z be nodal -- IK-ordered
%                                          y can be given as nodal or not.
%                                          y=f(x,z)
%                                          If y is constant, it is a 2D
%                                          surface. If y is not constant,
%                                          it can be interpreted as a 3D
%                                          surface.
%
% for surface in (y,z) plane, e.g. P(y,z), need y, z be nodal -- JK-ordered
%                                          x can be given as nodal or not.
%                                          x=f(y,z)
%                                          If x is constant, it can be
%                                          interpreted as a 2D surface. If
%                                          x is not constant, it can be
%                                          interpreted as a 3D surface
%
% for surface in 3D without order, this is not doable here. It can 
% be done with FEsurfaces, in which we do not need order
% 


   %
   %give location for each variable in this zone
   %
   NV=tdata.Nvar;
   for iv=1:NV   %NV is number of variables    
       if(tdata.surfaces(isurf).varloc~=1)
           var_location = 0;  %nodal 
       else    %variable iv is cell centered
           var_location = 1;  %cell centered
       end
       %make coordinate variables nodal 
       %assuming first variable is x, second is y, third is z
       %
       switch(surface_order)
           case 3                 %IJ-ordered surface data
               if(iv==1 || iv ==2)  %x,y must be nodal
                   var_location=0;
               end
           case 2                 %IK-ordered surface data
               if(iv==1 || iv ==3) %x,z must be nodal
                   var_location=0;
               end
           case 1                 %JK-ordered surface data
                                  %y,z  must be nodal
               if(iv==2 || iv ==3)
                   var_location=0;
               end
       end
       fwrite(fid_out,var_location,'int32');
   end
 end

%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%
%
%               
%
%
% Here we set it to false for ordered data (ordered surface)
% flase = 0 (int32), true = 1 (int32)
%
dummy_int32 = 0; %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,dummy_int32,'int32');

% For face connection (for each face in curren tzone), Tecplot takes the following 
% *****The combination of cell and face numbers in the current zone must be unique*****
% ***** multiple entries are not allowed                                          *****
%-----------------------------------------------------------------------------------------------------
% FaceNeighbor Mode |  #Values  | Data                                                               
%-----------------------------------------------------------------------------------------------------
% LocalOneToOne     |   3       |  cz,fz,cz   (cz--cell number in current zone                       
%                   |           |              fz--cell face number in current zone
%                   |           |              cz--cell number of the neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% LocalOneToMany    |   nz+4    |  cz,fz,oz,nz,cz1,cz2,...,czn
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--number of cell face in current zone
%                   |           |        oz--face obscuration flag
%                   |           |        nz--number of neighboring cells for one-to-many options
%                   |           |        cz1--cell number of the neighbor cell in current zone
%                   |           |        cz2--cell number of 2nd neighbor cell in current zone
%                   |           |             ...
%                   |           |        czn--cell number of the n'th neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToOne    |   4       |  cz,fz,ZZ,CZ
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--face number in current zone
%                   |           |        ZZ--remote zone number
%                   |           |        CZ--cell number of neighboring cell in remote zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToMany   |   2*nz+4  |  cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%-----------------------------------------------------------------------------------------------------
%
% cz --cell in current zone
% fz --face number of cell in current zone
% oz --face obscuration flag (only applies to one-to-many)
%         0 -- face partially obscured
%         1 -- face entirely obscured
% nz --number of cell or zone/cell associateions
% ZZ --remote Zone
% CZ --cell in remote zone
%
% cz,fz combinations must be unique. Additionally,Tecplot assume sthat with
% the one-to-one face neighbor modes a supplied cell face is entirely
% obscured by its neighbor. With one-to-many, the obscuration flag must be 
% supplied. 
%
% See doc/adkum.pdf page 71 for Face Neihbors 

%Use face neighbors to specify connections between zones (global connections) or connections within zones
%(local connections) for ordered or classic finite-element data1. Face neighbor connections are used when
%deriving variables or drawing contour lines. Specifying face neighbors, typically leads to smoother
%connections. NOTE: face neighbors have expensive performance implications. Use face neighbors to
%manually specify connections that are not defined via the connectivity list.

%The nature of the data arranged in the Face Neighbor Connections list depends upon the Face Neighbor
%Mode, described in the above table. To connect the cells along one edge to cells on another edge of the
%same zone, use one of the "local" types. To connect cells of one zone to cells of another zone or zones, use
%one of the "global" types. If the points of the cells are exactly aligned with the neighboring cell points, use
%"one-to-one". If even one cell face is neighbor to two or more other cell faces, use "one-to-many".
%


% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

%No misc user defined face neighbor connections for ordered surface data

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

if (NoOfUserDefinedNeighbourConn ~=0)  
    
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)

     %Pls modify based on data
     user_defined_face_neighbor_mode= 0; %Local 1-to-1
     
     
     fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
     
     %skiped below for ordered data 
     
%      if(zone_type>0)  %for non-oredered zones only
% %           if FE Zone:
% %             +-----------+
% %             | INT32     |     Indicates if the finite element face neighbors are
% %             +-----------+     completely specified by the miscellaneous face neighbors
% %                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
% %                               auto assignment of face neighbors otherwise all faces not
% %                               specified are considered boundaries. If no, then Tecplot will
% %                               perform auto-assignment of the face neighbors unless the
% %                               raw face neighbor array was supplied. This option is not
% %                               valid for ORDERED zones.
%         
%         %pls modify this value based on your data
%         fe_face_complete_by_misc_face_neighbors=1;
%         
%         fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
%      end
    
    
end

if(zone_type==0)  %for ordered zone

%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   % ordered zone, 3 int32s for num points i -j -k
   %
   
   fwrite(fid_out,Imax,'int32');
   fwrite(fid_out,Jmax,'int32');
   fwrite(fid_out,Kmax,'int32');

else  %for non-ordered (FE) zones
   
% %         if FE Zone:
% %           +-----------+
% %           | INT32     |     NumPts
% %           +-----------+
% % 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
% % 		        +-----------+ 
% % 			    | INT32     |   NumFaces 
% % 	    		+-----------+ 
% % 			    +-----------+ 
% % 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% % 			    +-----------+   zones, this is NumFaces*2. 
% % 			    +-----------+ 
% % 			    | INT32	    |   Total number of boundary faces. If any 
% % 			    +-----------+   boundary faces exist, include one to represent 
% % 					            no neighboring element. 
% % 			    +-----------+ 
% % 			    | INT32	    |   Total number of boundary connections. 
% % 			    +-----------+
% %
% %           +-----------+
% %           | INT32     |     NumElements.
% %           +-----------+       
% %           +-----------+
% %           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
% %           +-----------+       
% 
%      % %uncomment and complete here
% 
%       %NumPts=  ; %give number of nodes
%       %fwrite(fid_out,NumPts,'int32');
%       
         if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
% 		        +-----------+ 
% 			| INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32     |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      

%          %uncomment the following and complete
%          %
%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');
         end      

%       %NumElements= ; %number of elements
%       %fwrite(fid_out,Elements,'int32');
%       
%       ICellDim=0;
%       JCellDim=0;
%       KCellDim=0;
%       fwrite(fid_out,ICellDim,'int32');
%       fwrite(fid_out,JCellDim,'int32');
%       fwrite(fid_out,KCellDim,'int32');
    
end


% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%

%check if have aux data for this surface

have_auxdata=(   ~isempty(find(strcmp(surffields,'auxname')==1))  ...
              && ~isempty(find(strcmp(surffields,'auxval')==1)));   

if(have_auxdata)
   aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
   aux_names=tdata.surfaces(isurf).auxname;  %cell array or string
   aux_values=tdata.surfaces(isurf).auxval; %cell array or string
   
   %
   %check if they are cell arrays (1D) or just a single string
   %

    if(iscellstr(aux_names) && iscellstr(aux_values))  %cells       
     
      N_aux=min(length(aux_names),length(aux_values));
      %or check if they are string (not array), just single string value
      %number of aux is obtained from aux_names
      %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values) > N_aux)
          %give warnging if they do not match in size
          %warning(['Warning: surface ' int2str(isurf) ' auxname and auxval have inconsistent sizes']);
          %warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
    else  %single string
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
    end
  
   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: surface ' int2str(isurf) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this line']);
       aux_data=0;
   end  
else 
   aux_data=0;  %no aux data
end

%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.
if(aux_data==0)
   fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
        
        fwrite(fid_out,aux_data,'int32');    %repeat for each aux variable
        
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
               aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value  
    end
end

end  %end of isurf loop

%----
%Loop through all cubes (3D ordered volume data) and treat each cube as a zone
%I-J-K ordered data, used to describe a 3D volume, yet with data
%points known to follow I-J-K orders, if one of the dimension sizes is one,
%it degrades to ordered surfaces. If two of the dimension sizes are one,
%then it degrades to ordered lines
%-----

have_cubes=~isempty(find(strcmp(tdatanames,'cubes')==1));  %check if have cubes
%get number of cubes in tdata
if(have_cubes)  %make sure have cubes
    if(isstruct(tdata.cubes))     %make sure tdata.cubes is structure array
                                  %or at least structure 
       Ncubes =length(tdata.cubes);
    else
       Ncubes = 0;
    end
else
    Ncubes=0;   
end

for icube=1:Ncubes

cubefields=fieldnames(tdata.cubes(icube));  %obtain field names
                                            %in cubes(icube)    
% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name 'data zone' and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

%check of zonename is available
have_zonename=~isempty(find(strcmp(cubefields,'zonename')==1));
if(have_zonename && isempty(tdata.cubes(icube).zonename))
   have_zonename=0;
end
if(~have_zonename)  %if no zone name, give one by default
                    %as 'Cube1','Cube2' etc, depending on index icube
   zone_name = ['Cube' int2str(icube)];  
else
   if(   ischar(tdata.cubes(icube).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.cubes(icube).zonename))) %and not an empty string
       zone_name=tdata.cubes(icube).zonename; 
   else
       %give warning and make up zone name
       warning(['cube number ' int2str(icube) 'zonename is NOT a string!']);
       zone_name=['Cube' int2str(icube)];  %make up names by default
       warning(['set to default value as ' zone_name]);
   end
end
                          
plt_write_string(fid_out, zone_name);

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

%always set to -1 for no parent zones
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

%check if have strandID
have_strandID=~isempty(find(strcmp(cubefields,'strandID')==1)); 
if(have_strandID && isempty(tdata.cubes(icube).strandID))
    have_strandID=0;
end

if(~have_strandID) %if not exist, then give zero by default
   strand_ID= -2;
else
    if(isinteger(tdata.cubes(icube).strandID))  %make sure it is an integer
        strand_ID=tdata.cubes(icube).strandID;
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.cubes(icube).strandID);
        catch 
            strand_ID=-2;
        end
        warning(['strandID of cube ' int2str(icube) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end

fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

%check if have solutiontime
have_solutiontime=~isempty(find(strcmp(cubefields,'solutiontime')==1));
if(have_solutiontime && isempty(tdata.cubes(icube).solutiontime))
    have_solutiontime=0;
end
if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(  isnumeric(tdata.cubes(icube).solutiontime) ...
       &&  isfinite(tdata.cubes(icube).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.cubes(icube).solutiontime;
    else
        warning(['solutiontime of cube ' int2str(icube) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');


%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

%this is default to -1 and to be determined by tecplot
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE,
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON,
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON
%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%
%

%For cubes (IJK ordered volume data) this is set to 0
zone_type = 0;  %ordered 
fwrite(fid_out,zone_type,'int32');

%
%for ordered volume data, we need to make sure x,y,z all exist and 
%have the same size (IMAX,JMAX,KMAX)
%IF any of IMAX, JMAX, KMAX =1, we have to also make sure that 
%the rest of variables have same size and non-centered
%

have_x=~isempty(find(strcmp(cubefields,'x')==1));
have_y=~isempty(find(strcmp(cubefields,'y')==1));
have_z=~isempty(find(strcmp(cubefields,'z')==1));
if(have_x && isempty(tdata.cubes(icube).x))
   have_x=0;
end
if(have_y && isempty(tdata.cubes(icube).y))
   have_y=0;
end
if(have_z && isempty(tdata.cubes(icube).z))
   have_z=0;
end

%
%Check variable size to make sure coordinate variables are nodal and have same size
%
if(have_x && have_y && have_z)  %must have x, y, z and their sizes must match
   Imax_x=size(tdata.cubes(icube).x,1);
   Jmax_x=size(tdata.cubes(icube).x,2);
   Kmax_x=size(tdata.cubes(icube).x,3);

   Imax_y=size(tdata.cubes(icube).y,1);
   Jmax_y=size(tdata.cubes(icube).y,2);
   Kmax_y=size(tdata.cubes(icube).y,3);
   
   Imax_z=size(tdata.cubes(icube).z,1);
   Jmax_z=size(tdata.cubes(icube).z,2);
   Kmax_z=size(tdata.cubes(icube).z,3);
 
   if(   Imax_x== Imax_y && Imax_x == Imax_z && Imax_y ==Imax_z  ...
      && Jmax_x== Jmax_y && Jmax_x == Jmax_z && Jmax_y ==Jmax_z  ...
      && Kmax_x== Kmax_y && Kmax_x == Kmax_z && Kmax_y ==Kmax_z)
          %good, sizes match
          
          Imax=Imax_x;
          Jmax=Jmax_x;
          Kmax=Kmax_x;
          
   else
       s=-1;
       display(['Error: cube ' int2str(icube) ' does not have matching x, y, z coordinate dimensions!! ']);
       return;
   end
else %give error
   s=-1;
   display(['Error: cube ' int2str(icube) ' must have 3D coordinates x, y, z !! ']);
   return;
end


%
%Data packing
% 0- Black
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

% check if there is datapacking 
% have_datapacking=~isempty(find(strcmp(cubefields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.cubes(icube).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.cubes(icube).datapacking) ...
%       &&  isfinite(tdata.cubes(icube).datapacking) ...
%       &&  (  tdata.cubes(icube).datapacking ==0 ...
%           || tdata.cubes(icube).datapacking ==1))
%         data_packing=tdata.cubes(icube).datapacking;
%    else
%        warning(['datapacking of cube ' int2str(icube) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end

data_packing=0;
% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------


%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  not specifying
%   1 ----  Specify
%  

%
%Wen Long: note that when data_packing=1, i.e. point (IsBlock=0), var
%location must be nodal!! That is Cell Centered variables can only occur
%when data_packing is block. Cell Centered variables CAN'T happen when
%data_packing method is point. For Nodal variables, they can occur for both
%point packing and block packing methods.
%

%check if have varloc
have_varloc=~isempty(find(strcmp(cubefields,'varloc')==1));
if(have_varloc && isempty(tdata.cubes(icube).varloc))
   have_varloc=0 ;
end
if(~have_varloc) %give default value as 0 (nodal)
   var_loc=0;
   var_specified=0;
else
   if(     isfinite(tdata.cubes(icube).varloc) ...
      &&  (  tdata.cubes(icube).varloc ==0 ...
          || tdata.cubes(icube).varloc ==1))
        var_loc=tdata.cubes(icube).varloc;
        var_specified=1;
        
        %Wen Long, make sure when var_specified==1
        %data_packing is 0 (block) 
        %rule out conflict conditions (data_packing=1 (point) and var_loc=
        %1 (cell-center) cannot co-exist)
        %That is to say when var_loc=1, data_packing must be zero (block)
        %
        if(data_packing==1 && var_loc==1)
            s=-1;
            display(['Error: datapacking of cube ' int2str(icube) ...
                     ' (1- point) conflicts with varloc (1-centered)']);
            return;
        end
        
        %check if any dimension size is one, if yes, then
        %make sure we do not have cell-centered data. Volume cells do not
        %exist if any dimension size degrades to one
        %
          
        if(Imax_x==1 ||Jmax_x==1||Kmax_x==1)
              if(var_loc==1)
                s=-1;
                display(['Error: one dimension of  cube ' int2str(icube) ...
                       ' is singleton and cannot use cell-centered data']);
                return;
              end
        end
   else
       warning(['var location of cube ' int2str(icube) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
       var_loc=0;
       var_specified=0;
       warning(['set default value as 0 (nodal)']);
   end
end

 %output var_specified
 fwrite(fid_out,var_specified,'int32');

 if(var_specified==1)
%
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)
   %
   %give location for each variable in this zone
   %
      
   NV=tdata.Nvar;
   for iv=1:NV   %NV is number of variables    
       if(tdata.cubes(icube).varloc~=1)
           var_location = 0;  %nodal 
       else    %variable iv is cell centered
           var_location = 1; %cell centered
       end
       %make coordinate variables are given nodal 
       %assuming first variable is x, second is y, third is z
       %essentially, x, y, z coordiante variables in this cube
       %must be provided as nodal values
       %rest of variables can be centered if data_packing method is block
       %
       if(iv==1||iv==2||iv==3)
          var_location=0;
       end
       fwrite(fid_out,var_location,'int32');
   end
 
 end
 
%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%

% Here we set it to false for ordered data 
% flase = 0 (int32), true = 1 (int32)
%
dummy_int32 = 0; %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,dummy_int32,'int32');

% For face connection (for each face in curren tzone), Tecplot takes the following
% *****The combination of cell and face numbers in the current zone must be unique*****
% ***** multiple entries are not allowed                                          *****
%-----------------------------------------------------------------------------------------------------
% FaceNeighbor Mode |  #Values  | Data
%-----------------------------------------------------------------------------------------------------
% LocalOneToOne     |   3       |  cz,fz,cz   (cz--cell number in current zone
%                   |           |              fz--cell face number in current zone
%                   |           |              cz--cell number of the neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% LocalOneToMany    |   nz+4    |  cz,fz,oz,nz,cz1,cz2,...,czn
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--number of cell face in current zone
%                   |           |        oz--face obscuration flag
%                   |           |        nz--number of neighboring cells for one-to-many options
%                   |           |        cz1--cell number of the neighbor cell in current zone
%                   |           |        cz2--cell number of 2nd neighbor cell in current zone
%                   |           |             ...
%                   |           |        czn--cell number of the n'th neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToOne    |   4       |  cz,fz,ZZ,CZ
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--face number in current zone
%                   |           |        ZZ--remote zone number
%                   |           |        CZ--cell number of neighboring cell in remote zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToMany   |   2*nz+4  |  cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%-----------------------------------------------------------------------------------------------------
% cz --cell in current zone
% fz --face of cell in current zone
% oz --face obscuration flag (only applies to one-to-many)
%         0 -- face partially obscured
%         1 -- face entirely obscured
% nz --number of cell or zone/cell associateions
% ZZ --remote Zone
% CZ --cell in remote zone
%
% cz,fz combinations must be unique. Additionally,Tecplot assume sthat with
% the one-to-one face neighbor modes a supplied cell face is entirely
% obscured by its neighbor. With one-to-many, the obscuration flag must be 
% supplied. 
%

% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

if (NoOfUserDefinedNeighbourConn ~=0)  
    
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)

     %Pls modify based on data
     user_defined_face_neighbor_mode= 0; %Local 1-to-1
     
     
     fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
     
     if(zone_type>0)  %for non-oredered zones only
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
        
        %pls modify this value based on your data
        fe_face_complete_by_misc_face_neighbors=1;
        
        fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
     end
    
    
end

if(zone_type==0)  %for ordered zone

%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   % ordered zone var 3 int32s for num points i -j -k
   %
   
   %Imax, Jmax, Kmax already calculated above
   fwrite(fid_out,Imax,'int32');
   fwrite(fid_out,Jmax,'int32');
   fwrite(fid_out,Kmax,'int32');

else  %for non-ordered (FE) zones
   
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
% 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary connections. 
% 			    +-----------+
%
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       

     % %uncomment and complete here

      %NumPts=  ; %give number of nodes
      %fwrite(fid_out,NumPts,'int32');
         if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      

%          %uncomment the following and complete
%          %
%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');
         end      
      
      %NumElements= ; %number of elements
      %fwrite(fid_out,NumElements,'int32');
      
      ICellDim=0;
      JCellDim=0;
      KCellDim=0;
      fwrite(fid_out,ICellDim,'int32');
      fwrite(fid_out,JCellDim,'int32');
      fwrite(fid_out,KCellDim,'int32');
    
end

% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%

%=============
% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%
% Here we only deal with texts for 
%       tdata.cubes(icube).auxname and tdata.cubes(icube). 
%

%check if we have aux data for this cube, true if both auxname and 
%auxval exist and are 'cell arrays or string' and neither of them is empty
have_auxdata=(   ~isempty(find(strcmp(cubefields,'auxname')==1))  ...
              && ~isempty(find(strcmp(cubefields,'auxval')==1)));   

if(have_auxdata)

    aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
    aux_names=tdata.cubes(icube).auxname;  %cell array or string
    aux_values=tdata.cubes(icube).auxval;   %cell array or string
   %
   %check if they are cell arrays (1D) or just a single string
   %
   if(iscellstr(aux_names) && iscellstr(aux_values))
      %or check if they are string (not array), just single string value
      N_aux=min(length(aux_names),length(aux_values));
                             %number of aux is obtained from aux_names
                             %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values)> N_aux)
          %give warnging if they do not match in size
          warning(['Warning: cube ' int2str(icube) ' auxname and auxval have inconsistent sizes']);
          warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
   else
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
   end
   
   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: cube ' int2str(icube) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this cube']);
       aux_data=0;
   end
   
else
   aux_data=0;  %no aux data
end

%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.

if(aux_data==0)
    fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
            
          fwrite(fid_out,aux_data,'int32');  %repeat for each aux variable
          
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
              aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value
          
    end
end
end  %end of icube loop

%-----
%
%Loop trhough all Finite Element lines (line segement elements (cells)) 
% and treat each FEline as a zone. Each line element is made of two 
% end points variables can be defined on the element (either located 
% at the center of the line element (mid-point of the segment, 
% cell-center) or % at the end points (nodal) )
%
% Line element can be
%        1D (given by x coordiantes of the two end points)
%     or 2D (given by (x,y) coordiantes of the two end points)
%     or 3D (given by (x,y,z) coordiantes of the two end points)
% here, x, y, z are all 1D arrays
%
% Variables defined on the element is a collection of 1D variables
% if nodal, the size of the variable is equal to number of nodes
% if cell-centered, the size of the variable is equal to number of elements
%
% Nodes are defined by nodal coordinates (x) or (x,y) or (x,y,z)
% Elements are defined using nodal connectivit matrix. For each 
% element, two node numbers are assigned to give the compisition of nodes
% for the element
%

have_FElines=~isempty(find(strcmp(tdatanames,'FElines')==1));  %check if have FElines
%get number of FElines in tdata
if(have_FElines)  %make sure have FElines
    if(isstruct(tdata.FElines))  %make sure tdata.FElines is structure array
                                  %or at least structure 
       NFElines =length(tdata.FElines);
    else
       NFElines =0;
    end
else
    NFElines=0;   
end

for iFEline=1:NFElines
    
    FElinefields=fieldnames(tdata.FElines(iFEline));  %find fields in
                                             %FElines(iFEline) structure
%     
% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name 'data zone' and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

%check if have zonename
have_zonename=~isempty(find(strcmp(FElinefields,'zonename')==1));
if(have_zonename && isempty(tdata.FElines(iFEline).zonename))
   have_zonename=0;
end
if(~have_zonename)  %if no zone name, give one by default
                    %as 'FELine1','FELine2' etc, depending 
                    %on index iFEline
   zone_name = ['Line' int2str(iFEline)];  
else
   if(   ischar(tdata.FElines(iFEline).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.FElines(iFEline).zonename))) %and not an empty string
       zone_name=tdata.FElines(iFEline).zonename; 
   else
       %give warning and make up zone name
       warning(['FElines number ' int2str(iFEline) 'zonename is NOT a string!']);
       zone_name=['FELine' int2str(iFEline)];  %make up names by default
       warning(['set to default value as ' zone_name]);
   end
end
                          
plt_write_string(fid_out, zone_name);

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

%Always set to -1 for having no parent zone
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

%check if have strandID

have_strandID=~isempty(find(strcmp(FElinefields,'strandID')==1)); 
if(have_strandID && isempty(tdata.FElines(iFEline).strandID))
    have_strandID=0;
end

if(~have_strandID) %if not exist, then give zero by default
   strand_ID=-2;
else
    if(isinteger(tdata.FElines(iFEline).strandID))  %make sure it is an integer
        strand_ID=tdata.FElines(iFEline).strandID;
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.FElines(iFEline).strandID);
        catch 
            strand_ID=-2;
        end
        warning(['strandID of FEline ' int2str(iFEline) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end

fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

%check if this line has solutiontime
have_solutiontime=~isempty(find(strcmp(FElinefields,'solutiontime')==1));
if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(   isnumeric(tdata.FElines(iFEline).solutiontime) ...
       &&  isfinite(tdata.FElines(iFEline).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.FElines(iFEline).solutiontime;
    else
        warning(['solutiontime of FEline ' int2str(iFEline) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');

%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

%Always let tecplot choose zone color
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE,
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON,
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON
%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%
%

zone_type = 1;  %for FElines, zone_type is set to FELINESEG
fwrite(fid_out,zone_type,'int32');

%
%Data packing
% 0- Block
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

%check if have datapacking 
% have_datapacking=~isempty(find(strcmp(FElinefields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.FElines(iFEline).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.FElines(iFEline).datapacking) ...
%       &&  isfinite(tdata.FElines(iFEline).datapacking) ...
%       &&  (  tdata.FElines(iFEline).datapacking ==0 ...
%           || tdata.FElines(iFEline).datapacking ==1))
%         data_packing=tdata.FElines(iFEline).datapacking;
%    else
%        warning(['datapacking of FEline ' int2str(iFEline) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end
data_packing=0;
% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------



%
%Check variable size to make sure coordinate variables are nodal and have
%same size. x must have, y, z are optional. If y, z exist, they must
%have same size as x and must be 1D, size of x must be at list two points
%to make sure it can have at least one line element (segment)
%
%
%for 1D line, order of line = 0
%
%        x must exist, and the functions defined on this 1D line
% can be y=f(x), if y exists, y can be nodal or cell-centered
%        z=z(x), if z exists, z can be nodal or cell-centered
%        v=v(x), if v exists, v can be nodal or cell-centered
%
%
%for 2D line, it can be (depending on order/orientation of line)
%
%  case 1: order of line = 3
%             z=f(x,y),v=f(x,y), coordiantes x,y must exist
%             z is default to zeros if not given. 
%             v is optional. v is default to zeros if Nvar>3  
%             z and v can be nodal or cell-centered
%             x, y must be nodal
%  case 2: order of line = 2
%             y=f(x,z),v=f(x,z), coordiantes x,z must exist
%             y is default to zeros if not given
%             v is optional. v is default to zeros if Nvar>3  
%             y and v can be nodal or cell-centered
%             x,z must be nodal
%  case 3: order of line = 1
%             x=f(y,z),v=f(y,z), coordiantes y,z must exist
%             x is default to zeros if not given
%             v is optional. v is default to zeros if Nvar>3  
%             x and v can be nodal or cell-centered
%             y, z must be nodal
%
%for 3D line, order of line = 4
%
%  x, y, z must all exist, v must exist and v is 
%             v=f(x,y,z)
%             v can be nodal or cell-centered
%             x, y, z must be nodal

%
%by default order of line is 4, that is it is a 3D line
%

%check if have line order

have_order=~isempty(find(strcmp(FElinefields,'order')==1));
if(~have_order) %default to 4 (3D line)
    warning(['FEline ' int2str(iFEline) ' does not have line order specified']);
    FEline_order=4;
    warning(['set to default as ' int2str(FEline_order) ' (3D line)']);
else
    FEline_order=tdata.FElines(iFEline).order;
end

if(FEline_order<0 || FEline_order>4)
   s=-1;
   display(['Error: FEline ' int2str(iFEline) ' order is incorrect (must be 0,1,2,3 or 4)']);
   return;
end

%
%Check x,y,z and probe how the lines are defined
%
have_x=~isempty(find(strcmp(FElinefields,'x')==1));
have_y=~isempty(find(strcmp(FElinefields,'y')==1));
have_z=~isempty(find(strcmp(FElinefields,'z')==1));
if(have_x && isempty(tdata.FElines(iFEline).x))
   have_x=0;
end
if(have_y && isempty(tdata.FElines(iFEline).y))
   have_y=0;
end
if(have_z && isempty(tdata.FElines(iFEline).z))
   have_z=0;
end


%check variable sizes based on order
switch(FEline_order)
    case 0  %1D
        %make sure x exist, and size >=2
        if(have_x) 
            if(~isvector(tdata.FElines(iFEline).x))  
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' x value must be 1D array']);
                return;
            else
                Imax_x=length(tdata.FElines(iFEline).x(:)); 
                if(Imax_x<2)
                    s=-1;
                    display(['Error: FEline ' int2str(iFEline) ' x coordiante must have at least 2 points']);
                    return;
                else
                    NumNodes=Imax_x; %get number of nodes
                end
            end
        else
            s=-1;
            display(['Error: FEline ' int2str(iFEline) ' x must exist and be 1D array']);
            return;
        end
                
    case 1  %2D on X, Y coordinates
        %make sure x, y exist, size match and >=2
        if(have_x && have_y)
           DIMnumber_x=length(size(tdata.FElines(iFEline).x));
           DIMnumber_y=length(size(tdata.FElines(iFEline).y));
           %if(DIMnumber_x==1&&DIMnumber_y==1)
           if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).y))
               Imax_x=length(tdata.FElines(iFEline).x(:)); 
               Imax_y=length(tdata.FElines(iFEline).y(:)); 
               if(Imax_x~=Imax_y || Imax_x<2 || Imax_y<=2)
                   s=-1;
                   display(['Error: FEline ' int2str(iFEline) ' x y coordiante must have at least 2 points and size must match']);
                   return;
               else
                   NumNodes=Imax_x; %get number of nodes
               end
           else
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' x, y must be 1D arrays']);
                return
           end
            
        else
            s=-1;
            display(['Error: FEline ' int2str(iFEline) ' x and y must exist and be 1D array']);
            return
        end
        
    case 2  %2D on X, Z coordinates
        %make sure x, z exist, size match and >=2
        if(have_x && have_z)
           DIMnumber_x=length(size(tdata.FElines(iFEline).x));
           DIMnumber_z=length(size(tdata.FElines(iFEline).z));
           %if(DIMnumber_x==1&&DIMnumber_z==1)
           if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).z))
               Imax_x=length(tdata.FElines(iFEline).x(:)); 
               Imax_z=length(tdata.FElines(iFEline).z(:)); 
               if(Imax_x~=Imax_z || Imax_x<2 || Imax_z<2)
                   s=-1;
                   displat(['Error: FEline ' int2str(iFEline) ' x z coordiante must have at least 2 points and size must match']);
                   return
               else
                   NumNodes=Imax_x; %get number of nodes
               end
           else
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' x, z must be 1D arrays']);
                return
           end
            
        else
            s=-1;
            display(['Error: FEline ' int2str(iFEline) ' x and z must exist and be 1D array']);
            return
        end
    case 3  %2D on Y, Z coordiantes
        %make sure y, z exist, size match and >=2
        if(have_z && have_y)
           DIMnumber_z=length(size(tdata.FElines(iFEline).z));
           DIMnumber_y=length(size(tdata.FElines(iFEline).y));
           %if(DIMnumber_z==1&&DIMnumber_y==1)
           if(isvector(tdata.FElines(iFEline).z) && isvector(tdata.FElines(iFEline).y))
               Imax_z=length(tdata.FElines(iFEline).z(:)); 
               Imax_y=length(tdata.FElines(iFEline).y(:)); 
               if(Imax_z~=Imax_y || Imax_z<2 || Imax_y<2)
                   s=-1;
                   display(['Error: FEline ' int2str(iFEline) ' y z coordiante must have at least 2 points and size must match']);
                   return
               else
                   NumNodes=Imax_z; %get number of nodes
               end
           else
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' y, z must be 1D arrays']);
                return
           end  
        else
            s=-1;
            display(['Error: FEline ' int2str(iFEline) ' y and z must exist and be 1D array']);
            return
        end
    case 4  %3D on X, Y, Z coordiantes
        %make sure x, y, z exist, size match and >=2
        if(have_x && have_y && have_z)
           DIMnumber_x=length(size(tdata.FElines(iFEline).x));
           DIMnumber_y=length(size(tdata.FElines(iFEline).y));
           DIMnumber_z=length(size(tdata.FElines(iFEline).z));
           %if(DIMnumber_x==1&&DIMnumber_y==1&&DIMnumber_z==1)
           if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).y) ...
                                                 && isvector(tdata.FElines(iFEline).z))
               Imax_x=length(tdata.FElines(iFEline).x(:)); 
               Imax_y=length(tdata.FElines(iFEline).y(:)); 
               Imax_z=length(tdata.FElines(iFEline).z(:));  
               if(  Imax_x~=Imax_y ||Imax_x~=Imax_z ||Imax_y~=Imax_z  ...
                   || Imax_x<2 || Imax_y< 2||Imax_z< 2 )
                   s=-1;
                   display(['Error: FEline ' int2str(iFEline) ' x y z coordiante must have at least 2 points and size must match']);
                   return
               else
                   NumNodes=Imax_x; %get number of nodes
               end
           else
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' x, y, z must be 1D arrays']);
                return
           end
        else
            s=-1;
            display(['Error: FEline ' int2str(iFEline) ' x, y, z must exist and be 1D array']);
            return
        end
end

%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  not specifying
%   1 ----  Specify
%  

%check if varloc is specified
have_varloc=~isempty(find(strcmp(FElinefields,'varloc')==1));
if(have_varloc && isempty(tdata.FElines(iFEline).varloc))
   have_varloc=0;
end
if(~have_varloc) %give default value as 0 (nodal)
   var_loc=0;
   var_specified=0;
else
   if(     isfinite(tdata.FElines(iFEline).varloc) ...
      &&  (  tdata.FElines(iFEline).varloc ==0 ...
          || tdata.FElines(iFEline).varloc ==1))
        var_loc=tdata.FElines(iFEline).varloc;
        var_specified=1;
        
        %Wen Long, make sure when var_specified==1
        %data_packing is 0 (block) 
        %rule out conflict conditions (data_packing=1 (point) and var_loc=
        %1 (cell-center) cannot co-exist)
        %That is to say when var_loc=1, data_packing must be zero (block)
        %
        if(data_packing==1 && var_loc==1)
            s=-1;
            dispplay(['Error: datapacking of FEline ' int2str(iFEline) ' (1 point) conflicts with varloc']);
            return
        end
   else
       warning(['var location of FEline ' int2str(iFEline) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
       var_loc=0;
       var_specified=0;
       warning(['set default value as 0 (nodal)']);
   end
end

fwrite(fid_out,var_specified,'int32');

 if(var_specified==1)
%
%           +-----------+
%           | INT32*NV  |  Variable Location (only specify if above is 1).  
%           +-----------+  0 = Node, 1 = Cell Centered (See note 5.)
%
   %
   %give location for each variable in this zone
   %
   
   NV=tdata.Nvar;
   for iv=1:NV   %NV is number of variables 
       if(tdata.FElines(iFEline).varloc~=1)
           var_location = 0;  %nodal 
       else    %variable iv is cell centered
           var_location = 1; %cell centered
       end
       
       %make coordinate variables nodal
       %assuming first variable is x, second is y, third is z
       switch(FEline_order)
           
           case 0 %1D line defined on x
               if(iv==1)
                   var_location=0; %x must be nodal
               end
           case 1 %2D line defined on (x,y)
               if(iv==1||iv==2)  
                   var_location=0; %x and y must be nodal
               end
               
           case 2 %2D line defined on (x,z)
               if(iv==1||iv==3) 
                   var_location=0; %x and z must be nodal
               end
               
           case 3 %2D line defined on (y,z)
               if(iv==2||iv==3)
                   var_location=0; %y and z must be nodal
               end
           case 4 %3D line defined on (x,y,z)
               if(iv==1||iv==2||iv==3)
                   var_location=0; %x,y,and z must be nodal
               end
       end
       fwrite(fid_out,var_location,'int32');
   end
 
 end

%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face.
%                             Neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%

%
% Here we set it to false for FELINESEG zone
% flase = 0 (int32), true = 1 (int32)
%
raw_1to1_face_supplied=0;  %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,raw_1to1_face_supplied,'int32');

% For face connection (for each face in curren tzone), Tecplot takes the following
% *****The combination of cell and face numbers in the current zone must be unique*****
% ***** multiple entries are not allowed                                          *****
%-----------------------------------------------------------------------------------------------------
% FaceNeighbor Mode |  #Values  | Data
%-----------------------------------------------------------------------------------------------------
% LocalOneToOne     |   3       |  cz,fz,cz   (cz--cell number in current zone
%                   |           |              fz--cell face number in current zone
%                   |           |              cz--cell number of the neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% LocalOneToMany    |   nz+4    |  cz,fz,oz,nz,cz1,cz2,...,czn
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--number of cell face in current zone
%                   |           |        oz--face obscuration flag
%                   |           |        nz--number of neighboring cells for one-to-many options
%                   |           |        cz1--cell number of the neighbor cell in current zone
%                   |           |        cz2--cell number of 2nd neighbor cell in current zone
%                   |           |             ...
%                   |           |        czn--cell number of the n'th neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToOne    |   4       |  cz,fz,ZZ,CZ
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--face number in current zone
%                   |           |        ZZ--remote zone number
%                   |           |        CZ--cell number of neighboring cell in remote zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToMany   |   2*nz+4  |  cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%-----------------------------------------------------------------------------------------------------
% cz --cell in current zone
% fz --face of cell in current zone
% oz --face obscuration flag (only applies to one-to-many)
%         0 -- face partially obscured
%         1 -- face entirely obscured
% nz --number of cell or zone/cell associateions
% ZZ --remote Zone
% CZ --cell in remote zone
%
% cz,fz combinations must be unique. Additionally,Tecplot assume sthat with
% the one-to-one face neighbor modes a supplied cell face is entirely
% obscured by its neighbor. With one-to-many, the obscuration flag must be 
% supplied. 
%

% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

if (NoOfUserDefinedNeighbourConn ~=0)  
    
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)

     %Pls modify based on data
     user_defined_face_neighbor_mode= 0; %Local 1-to-1
     
     
     fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
     
     if(zone_type>0)  %for non-oredered zones only
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
        
        %pls modify this value based on your data
        fe_face_complete_by_misc_face_neighbors=1;
        
        fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
     end
    
    
end

if(zone_type==0)  %for ordered zone

%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   % ordered zone var 3 int32s for num points i -j -k
   %
   
   %Imax=
   %Jmax=
   %Kmax=
   
   %fwrite(fid_out,Imax,'int32');
   %fwrite(fid_out,Jmax,'int32');
   %fwrite(fid_out,Kmax,'int32');

else  %for non-ordered (FE) zones
   
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
% 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON:
% 
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	|   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	|   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					    no neighboring element. 
% 			    +-----------+ 
% 			    | INT32	|   Total number of boundary connections. 
% 			    +-----------+
% 
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       

     % %uncomment and complete here

     %find number of nodes and number of elements
      NumPts=NumNodes;  %NumNodes is calculated already
      fwrite(fid_out,NumPts,'int32');     

     %find NumFaces, total number of face nodes and total number
     %of boundary faces, total number of boundary connections
     %for FEPOLYGON and FEPOLYHEDRON
     
      if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. (for FEPOLYGON
%                               eacn face is a line segment made of two 
%                               nodes) (Wen Long: can this polygon be
%                               3D???, i.e. can polygon not belong
%                               to a plane surface?)
%
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      

%          %uncomment the following and complete
%          %
%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');
      end
      
      %find number of elements from e2n array of the line
      have_e2n=~isempty(find(strcmp(FElinefields,'e2n')==1));
      if(have_e2n && isempty(tdata.FElines(iFEline).e2n))
         have_e2n=0; 
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix have to be provided']);
          return;
      else
         NumElements=size(tdata.FElines(iFEline).e2n,1); %
         NumNodes_per_Element=size(tdata.FElines(iFEline).e2n,2); %
         %check if number of nodes per element is two
         if(NumNodes_per_Element~=2) %make sure two nodes per line
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix must have two columns']);
                return;
         end
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix must have at least 1 row (element)']);
             return;
         end
      end
      
      fwrite(fid_out,NumElements,'int32');
      
      %
      % Reserved by tecplot for future versions
      % ICellDim,JCellDim,KCellDim (for future use; set to zero)
      %
      
      ICellDim=0;
      JCellDim=0;
      KCellDim=0;
      fwrite(fid_out,ICellDim,'int32');
      fwrite(fid_out,JCellDim,'int32');
      fwrite(fid_out,KCellDim,'int32');
    
end


% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well

%check if we have aux data for this FEline, true if both auxname and 
%auxval exist and are 'cell arrays or string' and neither of them is empty
have_auxdata=(  ~isempty(find(strcmp(FElinefields,'auxname')==1))  ...
              && ~isempty(find(strcmp(FElinefields,'auxval')==1)));   
if(have_auxdata)
   aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
    aux_names=tdata.FElines(iFEline).auxname;  %cell array or string
   aux_values=tdata.FElines(iFEline).auxval;   %cell array or string
   %
   %check if they are cell arrays (1D) or just a single string
   %
   if(iscellstr(aux_names) && iscellstr(aux_values))
      %or check if they are string (not array), just single string value
      N_aux=min(length(aux_names),length(aux_values));
                             %number of aux is obtained from aux_names
                             %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values)> N_aux)
          %give warnging if they do not match in size
          warning(['Warning: FEline ' int2str(iFEline) ' auxname and auxval have inconsistent sizes']);
          warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
   else
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
   end
   
   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: FEline ' int2str(iFEline) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this FEline']);
       aux_data=0;
   end
   
else
   aux_data=0;  %no aux data
end

%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.
if(aux_data==0)
    fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
            
          fwrite(fid_out,aux_data,'int32');  %repeat for each aux variable
          
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
              aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value
          
    end
end

end

%-----------

%
%Loop through all FEsurfaces and treat each surface as a zone
%

have_FEsurfaces=~isempty(find(strcmp(tdatanames,'FEsurfaces')==1));  %check if have FEsurfaces
%get number of FEsurfaces in tdata
if(have_FEsurfaces)  %make sure have FEsurfaces
    if(isstruct(tdata.FEsurfaces))  %make sure tdata.surfaces is structure array
                                  %or at least structure 
       NFEsurfs =length(tdata.FEsurfaces);
    else
       NFEsurfs = 0;
    end
else
    NFEsurfs=0;   
end


for iFEsurf=1:NFEsurfs

    FEsurffields=fieldnames(tdata.FEsurfaces(iFEsurf));  %find fields in
                                                         %this FEsurface
%        +-----------+ 
%        | FLOAT32   |    Zone marker. Value = 299.0 
%        +-----------+ 
% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name 'data zone' and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

%check if have zonename
have_zonename=~isempty(find(strcmp(FEsurffields,'zonename')==1));
if(have_zonename && isempty(tdata.FEsurfaces(iFEsurf).zonename))
   have_zonename=0;
end
if(~have_zonename)  %if no zone name, give one by default
                    %as 'FESurface1','FESurface2' etc, depending 
                    %on index iFEsurf
   zone_name = ['FESurface' int2str(iFEsurf)];  
else
   if(      ischar(tdata.FEsurfaces(iFEsurf).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.FEsurfaces(iFEsurf).zonename))) %and not an empty string
       zone_name=tdata.FEsurfaces(iFEsurf).zonename; 
   else
       %give warning and make up zone name
       warning(['FEsurfaces number ' int2str(iFEsurf) 'zonename is NOT a string!']);
       zone_name=['FEsurface' int2str(iFEsurf)];  %make up names by default
       warning(['set to default value as ' zone_name]);
   end
end
                          
plt_write_string(fid_out, zone_name);

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

%Always set to -1 for having no parent zone
dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');


% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

have_strandID=~isempty(find(strcmp(FEsurffields,'strandID')==1)); 
if(have_strandID && isempty(tdata.FEsurfaces(iFEsurf).strandID))
    have_strandID=0;
end
if(~have_strandID) %if not exist, then give zero by default
   strand_ID= -2;
else
    if(isinteger(tdata.FEsurfaces(iFEsurf).strandID))  %make sure it is an integer
        strand_ID=tdata.FEsurfaces(iFEsurf).strandID;
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.FEsurfaces(iFEsurf).strandID);
        catch 
            strand_ID=-2;
        end
        warning(['strandID of FEsurface ' int2str(iFEsurf) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end

fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

%check if this FEsurface has solutiontime
have_solutiontime=~isempty(find(strcmp(FEsurffields,'solutiontime')==1));
if(have_solutiontime && isempty(tdata.FEsurfaces(iFEsurf).solutiontime))
    have_solutiontime=0;
end
if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(   isnumeric(tdata.FEsurfaces(iFEsurf).solutiontime) ...
       &&  isfinite(tdata.FEsurfaces(iFEsurf).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.FEsurfaces(iFEsurf).solutiontime;
    else
        warning(['solutiontime of FEsurface ' int2str(iFEsurf) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');


%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON

%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%
%

%For FEsurface, we can have triangle (2), quadrilateral (3), and polygon
%(6)

%find out the zone_type for FEsurface. 

      %find number of elements from e2n array of the surface
      have_e2n=~isempty(find(strcmp(FEsurffields,'e2n')==1));
      if(have_e2n && isempty(tdata.FEsurfaces(iFEsurf).e2n))
         have_e2n=0; 
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                   'connectivity) matrix/structure have to be provided']);
          return;
          
      else
                  
         %determine if e2n is a structure or array
         e2n_is_struct=isstruct(tdata.FEsurfaces(iFEsurf).e2n);
         
         if(~e2n_is_struct)  %not a struct
         
            %if array, then find about its number of columns, if 3--triangle,
            %then it is traingular, if 4 then it is quadrilateral, if else
            %give error
            
            DIMnumber=length(size(tdata.FEsurfaces(iFEsurf).e2n));
            if(DIMnumber~=2) %make sure e2n is 2D array
                s=-1;
                display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                         'connectivity) matrix/structure have to be 2D array']); 
                return;
            else  

               %get number of elements
               NumElements=size(tdata.FEsurfaces(iFEsurf).e2n,1);

               %get number of nodes per element
               NumNodes_per_element=size(tdata.FEsurfaces(iFEsurf).e2n,2);

               switch(NumNodes_per_element)
                   case 3
                       zone_type=2; %    2=FETRIANGLE
                   case 4
                       zone_type=3; %    3=FEQUADRILATERAL
                   otherwise
                       s=-1;
                       display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                         'connectivity) array must have 3 (triangle) or 4 (quadrilateral) columns']); 
                       return;
               end
               
            end
         
         else %if is a struct then check size of e2n.nodes to make sure
              %it is at least 3.
                            
               if(isvector(tdata.FEsurfaces(iFEsurf).e2n)) %make sure
                                                           %it is a
                                                           %structure
                                                           %vector
                   %get number of elements
                   NumElements=length(tdata.FEsurfaces(iFEsurf).e2n);
                   
                   %check all elements and make sure each element has at
                   %least 3 nodes
                   for icell=1:NumElements
                       %find field names in e2n. Must have nodes field
                       e2nfields=fieldnames(tdata.FEsurfaces(iFEsurf).e2n); 
                       have_nodes=~isempty(find(strcmp(e2nfields,'nodes')==1));
                       if(have_nodes && isempty(tdata.FEsurfaces(iFEsurf).e2n.node))
                          have_nodes=0;
                       end
                       if(have_nodes)                           
                           %make sure nodes is 1D vector array and also
                           %number of them is greater than 3 (to be
                           %polygon)
                           if(    isvector(tdata.FEsurfaces(iFEsurf).e2n.nodes) ...
                               &&   length(tdata.FEsurfaces(iFEsurf).e2n.nodes) >=3)
                               %do nothing
                           else
                               s=-1;
                               display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                                    'connectivity) nodes must be 1D array of at least 3 node numbers ']); 
                               return;
                           end
                       
                       else %without nodes 
                           s=-1;
                           display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                                    'connectivity) have no valid nodes ']); 
                           return;
                           
                       end
                   end
                   
                   zone_type=6;  %FEPOLYGON if passed all above checks
                                 %in this case, because a polygon does have fixed number of vertexes
                                 %and edges, user have to also define the faces of the polygon (edges)
                                 %for each of the element. Hence there will have to be a 
                                 %face2element mapping
                                 %For simple cases where a face is only neighboring one elemment, 
                                 %this is easier. For cases where a face is neighboring more then one element
                                 %the faceneighboring relationship has to also be provided.
 
                                 %When the faceneighbor element is from a different zone, then it 
                                 %also needs information of the zone that the face is neighboring to.
                                 %this is not yet implemented. 


                   
               else
                   s=-1;
                   display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                     'connectivity) must be a structure vector']); 
                   return;
               end
         end
          
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                      'connectivity) matrix/struct must have at least 1 (element)']);
             return;
         end
      end

fwrite(fid_out,zone_type,'int32');


%
%Data packing
% 0- Block
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

%check if have datapacking
% have_datapacking=~isempty(find(strcmp(FEsurffields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.FEsurfaces(iFEsurf).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.FEsurfaces(iFEsurf).datapacking) ...
%       &&   isfinite(tdata.FEsurfaces(iFEsurf).datapacking) ...
%       &&  (  tdata.FEsurfaces(iFEsurf).datapacking ==0 ...
%           || tdata.FEsurfaces(iFEsurf).datapacking ==1))
%         data_packing=tdata.FEsurfaces(iFEsurf).datapacking;
%    else
%        warning(['datapacking of FEsurface ' int2str(iFEsurf) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end
data_packing=0;

% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------


%find out the order (orientation/placement) of the FE surface
%and number of nodes 

have_order=~isempty(find(strcmp(FEsurffields,'order')==1));
if(have_order && isempty(tdata.FEsurfaces(iFEsurf).order))
   have_order=0; 
end
if(~have_order) %default to 4 (3D surface)
    warning(['FEsurface ' int2str(iFEsurf) ' does not have surface order specified']);
    FEsurf_order=4;
    warning(['set to default as ' int2str(FEsurf_order) ' (3D surface)']);
else
    FEsurf_order=tdata.FEsurfaces(iFEsurf).order;
end

if(FEsurf_order<1 || FEsurf_order>4)
   s=-1;
   display(['Error: FEsurface ' int2str(iFEsurf) ' order is incorrect (must be 1,2,3 or 4)']);
   return;
end

%
%check avaibility of x,y,z (coordinate variables) and size consistency
%based on order of surface. Also calculate number of nodes based on
%coordinate variables
%

have_x=~isempty(find(strcmp(FEsurffields,'x')==1));
have_y=~isempty(find(strcmp(FEsurffields,'y')==1));
have_z=~isempty(find(strcmp(FEsurffields,'z')==1));
if(have_x && isempty(tdata.FEsurfaces(iFEsurf).x))
   have_x=0;
end
if(have_y && isempty(tdata.FEsurfaces(iFEsurf).y))
   have_y=0;
end
if(have_z && isempty(tdata.FEsurfaces(iFEsurf).z))
   have_z=0;
end


switch(FEsurf_order)
    
    case 3  %surface is defined on (x,y), e.g. z=f(x,y), v=v(x,y)
            %if z does not exist, it is set to zero
            %(x,y) are coordinate variables, must be nodal
            
            if(have_x && have_y)
                
                %make sure x and y are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                          length(tdata.FEsurfaces(iFEsurf).y))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' x and y must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x and y provided on all nodes ']);
               return;
            end
        
    case 2  %surface is defined on (x,z), e.g. y=f(x,z),v=v(x,z)
            %if y does not exist, it is set to zero
            %(x,z) are coordinate variables, must be nodal
            if(have_x && have_z)
                
                %make sure x and z are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                          length(tdata.FEsurfaces(iFEsurf).z))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' x and z must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x and z provided on all nodes ']);
               return;
            end
            
    case 1  %surface is defined on (y,z), e.g. x=f(y,z),v=v(y,z)
            %if x does not exist, it is set to zero
            %(y,z) are coordinate variables, must be nodal
            
             if(have_z && have_y)
                
                %make sure x and y are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && length(tdata.FEsurfaces(iFEsurf).z)== ...
                          length(tdata.FEsurfaces(iFEsurf).y))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).z);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' y and z must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1; 
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' y and z provided on all nodes ']);
               return;
            end

    case 4  %surface is defined on 3D curving surface (x,y,z)
            %v=v(x,y,z)
            %(x,y,z) are coordinate vairables, must be nodal
            
            if(have_x && have_y && have_z)
                
                %make sure x, y and z are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                           length(tdata.FEsurfaces(iFEsurf).y) ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                           length(tdata.FEsurfaces(iFEsurf).z))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' x, y and z must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x, y and z provided on all nodes ']);
               return;
            end

end

%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  not specifying
%   1 ----  Specify
%  


%check if varloc is specified
have_varloc=~isempty(find(strcmp(FEsurffields,'varloc')==1));
if(have_varloc && isempty(tdata.FEsurfaces(iFEsurf).varloc))
   have_varloc=0;
end
if(~have_varloc) %give default value as 0 (nodal)
   var_loc=0;
   var_specified=0;
else
   if(     isfinite(tdata.FEsurfaces(iFEsurf).varloc) ...
      &&  (  tdata.FEsurfaces(iFEsurf).varloc ==0 ...
          || tdata.FEsurfaces(iFEsurf).varloc ==1))

        var_loc=tdata.FEsurfaces(iFEsurf).varloc;
        var_specified=1;
        
        %Wen Long, make sure when var_specified==1
        %data_packing is 0 (block) 
        %rule out conflict conditions (data_packing=1 (point) and var_loc=
        %1 (cell-center) cannot co-exist)
        %That is to say when var_loc=1, data_packing must be zero (block)
        %
        if(data_packing==1 && var_loc==1)
            s=-1;
            dispplay(['Error: datapacking of FEsurface ' int2str(iFEsurf) ' (1 point) conflicts with varloc']);
            return
        end
   else
       warning(['var location of FEsurface ' int2str(iFEsurf) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
       var_loc=0;
       var_specified=0;
       warning(['set default value as 0 (nodal)']);
   end
end

 fwrite(fid_out,var_specified,'int32');

 if(var_specified==1)
%
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)

   %
   %give location for each variable in this zone
   %
   
   NV=tdata.Nvar;
   for iv=1:NV   %NV is number of variables 
       if(tdata.FEsurfaces(iFEsurf).varloc~=1)
           var_location = 0;  %nodal 
       else    %variable iv is cell centered
           var_location = 1; %cell centered
       end
       
       %make coordinate variables are nodal
       %assuming first variable is x, second is y, third is z
       %
       switch(FEsurf_order)
                     
           case 3 %2D surface defined on (x,y)
               if(iv==1||iv==2)  
                   var_location=0; %x and y must be nodal
               end
               
           case 2 %2D surface defined on (x,z)
               if(iv==1||iv==3) 
                   var_location=0; %x and z must be nodal
               end
               
           case 1 %2D surface defined on (y,z)
               if(iv==2||iv==3)
                   var_location=0; %y and z must be nodal
               end
           case 4 %3D surface defined on (x,y,z)
               if(iv==1||iv==2||iv==3)
                   var_location=0; %x,y,and z must be nodal
               end
       end
       fwrite(fid_out,var_location,'int32');
   end
   
 end
 
%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%

% Here we set it to false for ordered data 
% flase = 0 (int32), true = 1 (int32)
%
 %raw local 1-to-1 face neighbors not supplied
raw_1to1_face_supplied=0;  %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,raw_1to1_face_supplied,'int32');

% For face connection (for each face in curren tzone), Tecplot takes the following
% *****The combination of cell and face numbers in the current zone must be unique*****
% ***** multiple entries are not allowed                                          *****
%-----------------------------------------------------------------------------------------------------
% FaceNeighbor Mode |  #Values  | Data
%-----------------------------------------------------------------------------------------------------
% LocalOneToOne     |   3       |  cz,fz,cz   (cz--cell number in current zone
%                   |           |              fz--cell face number in current zone
%                   |           |              cz--cell number of the neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% LocalOneToMany    |   nz+4    |  cz,fz,oz,nz,cz1,cz2,...,czn
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--number of cell face in current zone
%                   |           |        oz--face obscuration flag
%                   |           |        nz--number of neighboring cells for one-to-many options
%                   |           |        cz1--cell number of the neighbor cell in current zone
%                   |           |        cz2--cell number of 2nd neighbor cell in current zone
%                   |           |             ...
%                   |           |        czn--cell number of the n'th neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToOne    |   4       |  cz,fz,ZZ,CZ
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--face number in current zone
%                   |           |        ZZ--remote zone number
%                   |           |        CZ--cell number of neighboring cell in remote zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToMany   |   2*nz+4  |  cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%-----------------------------------------------------------------------------------------------------

%
% cz --cell in current zone
% fz --face of cell in current zone
% oz --face obscuration flag (only applies to one-to-many)
%         0 -- face partially obscured
%         1 -- face entirely obscured
% nz --number of cell or zone/cell associateions
% ZZ --remote Zone
% CZ --cell in remote zone
%
% cz,fz combinations must be unique. Additionally,Tecplot assume sthat with
% the one-to-one face neighbor modes a supplied cell face is entirely
% obscured by its neighbor. With one-to-many, the obscuration flag must be 
% supplied. 
%

% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

if (NoOfUserDefinedNeighbourConn ~=0)  
    
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)

     %Pls modify based on data
     user_defined_face_neighbor_mode= 0; %Local 1-to-1
     
     
     fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
     
     if(zone_type>0)  %for non-oredered zones only
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
        
        %pls modify this value based on your data
        fe_face_complete_by_misc_face_neighbors=1;
        
        fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
     end
end

if(zone_type==0)  %for ordered zone

%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   % ordered zone var 3 int32s for num points i -j -k
   %
   
%    Imax=
%    Jmax=
%    Kmax=
%    fwrite(fid_out,Imax,'int32');
%    fwrite(fid_out,Jmax,'int32');
%    fwrite(fid_out,Kmax,'int32');

else  %for non-ordered (FE) zones
   
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
%                     if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
%                           +-----------+ 
%                           | INT32     |   NumFaces 
%                           +-----------+ 
%                           +-----------+ 
%                           | INT32     |   Total number of face nodes. For FEPOLYGON 
%                           +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	|   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					    no neighboring element. 
% 			    +-----------+ 
% 			    | INT32	|   Total number of boundary connections. 
% 			    +-----------+
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       

     % %uncomment and complete here

      NumPts=NumNodes  ; %give number of nodes
      fwrite(fid_out,NumPts,'int32');

         if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
%                             +-----------+ 
%                             | INT32     |   NumFaces 
%                             +-----------+ 
%                             +-----------+ 
%                             | INT32     |   Total number of face nodes. For FEPOLYGON 
%                             +-----------+   zones, this is NumFaces*2. 
%                             +-----------+ 
%                             | INT32     |   Total number of boundary faces. If any 
%                             +-----------+   boundary faces exist, include one to represent 
%                                                     no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      
%
%          %uncomment the following and complete
%          %

              % %(This is not yet implemented)
              % %NumFaces, NumFaceNodes, NumBryFaces, NumBryConnections
              % %etc information can be provided through
              % %e2n.Faces e2n.BounaryFace etc  
              % %

%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');
         end      
      
      
      %NumElements already calculated above
      
      fwrite(fid_out,NumElements,'int32');
      
      %Reserved by tecplot for future use; set to zero here
      ICellDim=0;
      JCellDim=0;
      KCellDim=0;
      fwrite(fid_out,ICellDim,'int32');
      fwrite(fid_out,JCellDim,'int32');
      fwrite(fid_out,KCellDim,'int32');
    
end


% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%

%check if we have aux data for this FEsurface, true if both auxname and 
%auxval exist and are 'cell arrays or string' and neither of them is empty
have_auxdata=(   ~isempty(find(strcmp(FEsurffields,'auxname')==1))  ...
              && ~isempty(find(strcmp(FEsurffields,'auxval')==1)));   
if(have_auxdata)
   aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
    aux_names=tdata.FEsurfaces(iFEsurf).auxname;  %cell array or string
   aux_values=tdata.FEsurfaces(iFEsurf).auxval;   %cell array or string
   %
   %check if they are cell arrays (1D) or just a single string
   %
   if(iscellstr(aux_names) && iscellstr(aux_values))
      %or check if they are string (not array), just single string value
      N_aux=min(length(aux_names),length(aux_values));
                             %number of aux is obtained from aux_names
                             %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values)> N_aux)
          %give warnging if they do not match in size
          warning(['Warning: FEsurface ' int2str(iFEsurf) ' auxname and auxval have inconsistent sizes']);
          warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
   else
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
   end
   
   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: FEsurface ' int2str(iFEsurf) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this FEsurface']);
       aux_data=0;
   end
   
else
   aux_data=0;  %no aux data
end

%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.
if(aux_data==0)
    fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
            
          fwrite(fid_out,aux_data,'int32');  %repeat for each aux variable
          
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
              aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value
          
    end
end
end

%------
%
% Loop through all FEvolumes and treat each volume as a zone
%
%--------
have_FEvolumes=~isempty(find(strcmp(tdatanames,'FEvolumes')==1));  %check if have FEvolume
%get number of FEvolumes in tdata
if(have_FEvolumes)  %make sure have FEvolume
    if(isstruct(tdata.FEvolumes))  %make sure tdata.surfaces is structure array
                                  %or at least structure 
       NFEvolumes =length(tdata.FEvolumes);
    else
       NFEvolumes = 0;
    end
else
    NFEvolumes =0;   
end

for iFEvol=1:NFEvolumes

    FEvolfields=fieldnames(tdata.FEvolumes(iFEvol));  %find fields in this 
                                                      %finite elment volume
                                                      %zone
% 
% write zone marker float32 = 299.0
% zone starts with a zone maker called 299.0
%

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% write zone name 'data zone' and null
%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+

have_zonename=~isempty(find(strcmp(FEvolfields,'zonename')==1));
if(have_zonename && isempty(tdata.FEvolumes(iFEvol).zonename))
   have_zonename=0;
end
if(~have_zonename)  %if no zone name, give one by default
                    %as 'FEVolume1','FEVolume2' etc, depending 
                    %on index iFEvol
   zone_name = ['FEVolume' int2str(iFEvol)];  
else
   if(      ischar(tdata.FEvolumes(iFEvol).zonename) ...   %make sure it is a string
      &&  (~isempty(tdata.FEvolumes(iFEvol).zonename))) %and not an empty string
       zone_name=tdata.FEvolumes(iFEvol).zonename; 
   else
       %give warning and make up zone name
       warning(['FEvolumes number ' int2str(iFEvol) 'zonename is NOT a string!']);
       zone_name=['FEvolume' int2str(iFEvol)];  %make up names by default
       warning(['set to default value as ' zone_name]);
   end
end
                          
plt_write_string(fid_out, zone_name);

% parent zone = 0 (no parent zone) int32
%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.

dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

% strand ID suggested as zero from old calls
%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID

%check if have strandID
have_strandID=~isempty(find(strcmp(FEvolfields,'strandID')==1)); 
if(have_strandID && isempty(tdata.FEvolumes(iFEvol).strandID))
    have_strandID=0;
end
if(~have_strandID) %if not exist, then give zero by default
   strand_ID= -2;
else
    if(isinteger(tdata.FEvolumes(iFEvol).strandID))  %make sure it is an integer
        strand_ID=tdata.FEvolumes(iFEvol).strandID;
    else
        try   %see if can conver to integer
            strand_ID=int32(tdata.FEvolumes(iFEvol).strandID);
        catch 
            strand_ID= -2;
        end
        warning(['strandID of FEvolume ' int2str(iFEvol) ' is not integer!!!']);
        warning(['set to default value as ' int2str(strand_ID)]);
    end
end

fwrite(fid_out,strand_ID,'int32');

%
% solution_time
%
%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
%

have_solutiontime=~isempty(find(strcmp(FEvolfields,'solutiontime')==1));
if(have_solutiontime && isempty(tdata.FEvolumes(iFEvol).solutiontime))
    have_solutiontime=0;
end
if(~have_solutiontime)
    solution_time=0; %defaut to zero if not given
else
    if(   isnumeric(tdata.FEvolumes(iFEvol).solutiontime) ...
       &&  isfinite(tdata.FEvolumes(iFEvol).solutiontime))
        %make sure it is numerica value (integer or float or double)
        %and it is finite
        solution_time=tdata.FEvolumes(iFEvol).solutiontime;
    else
        warning(['solutiontime of FEvolume ' int2str(iFEvol) ' is not numeric!!!']);
        solution_time=0;
        warning(['set to default value as ' num2str(solution_time)]);
    end
end
fwrite(fid_out,solution_time,'float64');

%Zone color 
% not used set to -1 (int32)
%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).

dummy_int32 = -1;
fwrite(fid_out,dummy_int32,'int32');

%
% zone_type 
%
%         +-----------+
%         | INT32     |   ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+   3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                         6=FEPOLYGON,7=FEPPOLYHEDRON
%
%  Teplot takes the following zone types:
%
%    0=ORDERED
%    1=FELINESEG
%    2=FETRIANGLE,
%    3=FEQUADRILATERAL
%    4=FETETRAHEDRON,
%    5=FEBRICK
%    6=FEPOLYGON
%    7=FEPOLYHEDRON
%
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%

%
%find zone_type based on e2n of this volume. If e2n is array of size NEx4
%then it is  4=FETETRAHEDRON. If e2n is array of size NEx8, then it is
%5=FEBRICK. Otherwise, it is 7=FEPOLYHEDRON and faces and element to node 
%connectivity must be provided in e2n as a structure
%

      %find number of elements from e2n array of the line
      have_e2n=~isempty(find(strcmp(FEvolfields,'e2n')==1));
      if(have_e2n && isempty(tdata.FEvolumes(iFEvol).e2n))
         have_e2n=0;
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                   'connectivity) matrix/structure have to be provided']);
          return;
          
      else
                  
         %determine if e2n is a structure or array
         e2n_is_struct=isstruct(tdata.FEvolumes(iFEvol).e2n);
         if(~e2n_is_struct)  %not a struct
         
            %if array, then find about its number of columns, if 3--triangle,
            %then it is traingular, if 4 then it is quadrilateral, if else
            %give error
            
            DIMnumber=length(size(tdata.FEvolumes(iFEvol).e2n));
            if(DIMnumber~=2) %make sure e2n is 2D array
                s=-1;
                display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                         'connectivity) matrix/structure have to be 2D array']); 
                return;
            else  

               %get number of elements
               NumElements=size(tdata.FEvolumes(iFEvol).e2n,1);

               %get number of nodes per element
               NumNodes_per_element=size(tdata.FEvolumes(iFEvol).e2n,2);

               switch(NumNodes_per_element)
                   case 4
                       zone_type=4; %    4=FETETRAHEDRON
                   case 8
                       zone_type=5; %    5=FEBRICK
                   otherwise
                       s=-1;
                       display(['Error: FEsurface ' int2str(iFEvol) ' e2n (element to node ' ...
                         'connectivity) array must have 4 (tetrahedron) or 8 (brick) columns']); 
                       return;
               end
               
            end
         
         else %if is a struct then check size of e2n.nodes to make sure
              %it is at least 4.
                            
               if(isvector(tdata.FEvolumes(iFEvol).e2n)) %make sure
                                                         %it is a
                                                         %structure
                                                         %vector
                   %get number of elements
                   NumElements=length(tdata.FEvolumes(iFEvol).e2n);
                   
                   %check all elements and make sure each element has at
                   %least 4 nodes 
                   
                   for icell=1:NumElements
                       %find field names in e2n. Must have nodes field
                       e2nfields=fieldnames(tdata.FEvolumes(iFEvol).e2n); 
                       have_nodes=~isempty(find(strcmp(e2nfields,'nodes')==1));
                       if(have_nodes && isempty(tdata.FEvolumes(iFEvol).e2n.nodes))
                          have_nodes=0;
                       end
                       if(have_nodes)
                           %make sure nodes is 1D vector array and also
                           %number of them is greater than 4 (to be
                           %tetrahedron)
                           
                           if(    isvector(tdata.FEvolumes(iFEvol).e2n.nodes)  ...
                               && length(tdata.FEvolumes(iFEvol).e2n.nodes) >=4)
                               %do nothing
                           else
                               s=-1;
                               display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                                    'connectivity) nodes must be 1D array of at least 4 node numbers ']); 
                               return;
                           end
                       
                       else %without nodes for this element
                           s=-1;
                           display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                                    'connectivity) have no valid nodes ']); 
                           return;
                           
                       end
                   end
                   
                   zone_type=7;  %FEPOLYHEDRON if passed all above checks

                                 %For FEPOLYHEDRON, each polyhedron element does not 
                                 %have to have fixed number of vertexes (nodes) and faces
                                 %it is important that the user also provide definition
                                 %of all faces of the polyhedron. This 
                                 %can be done with a face2node mapping matrix
                                 %listing all faces and their consisting nodes. The list
                                 %of each face must be in all clock-wise order or counter-clockwise order
 
                                 %Further more, all faces can be neighboring one or more than 
                                 %one element of the current zone or other zones 
                                 %the face neighboring mapping has to be also provided.

               else
                   s=-1;
                   display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                     'connectivity) must be a structure vector']); 
                   return;
               end
         end
          
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                   'connectivity) matrix/struct must have at least 1 (element)']);
             return;
         end
      end

fwrite(fid_out,zone_type,'int32');

%
%Data packing
% 0- Block
% 1- Point
%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+               

%check if have datapacking
% have_datapacking=~isempty(find(strcmp(FEvolfields,'datapacking')==1));
% if(have_datapacking && isempty(tdata.FEvolumes(iFEvol).datapacking))
%    have_datapacking=0;
% end
% if(~have_datapacking) %give default value as block (0)
%    data_packing=0;
% else
%    if(   isinteger(tdata.FEvolumes(iFEvol).datapacking) ...
%       &&   isfinite(tdata.FEvolumes(iFEvol).datapacking) ...
%       &&  (  tdata.FEvolumes(iFEvol).datapacking ==0 ...
%           || tdata.FEvolumes(iFEvol).datapacking ==1))
%         data_packing=tdata.FEvolumes(iFEvol).datapacking;
%    else
%        warning(['datapacking of FEvolume ' int2str(iFEvol) ' is neither 0 (block) nor 1 (point)!!!']);
%        data_packing=0;
%        warning(['set default value as 0 (block)']);
%    end
% end

data_packing=0;
% %output data_packing
% ---Wen Long- This is now deprecated and deleted----
% fwrite(fid_out,data_packing,'int32');
% ------------------------------------------------


%
% Whether or not specify variable location
%

% +-----------+
% | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
% +-----------+       located at the nodes.  1 = Specify
%
%   0 ----  not specifying
%   1 ----  Specify
%  


%check if varloc is specified
have_varloc=~isempty(find(strcmp(FEvolfields,'varloc')==1));
if(have_varloc && isempty(tdata.FEvolumes(iFEvol).varloc))
   have_varloc=0;
end
if(~have_varloc) %give default value as 0 (nodal)
   var_loc=0;
   var_specified=0;
else
   if(     isfinite(tdata.FEvolumes(iFEvol).varloc) ...
      &&  (  tdata.FEvolumes(iFEvol).varloc ==0 ...
          || tdata.FEvolumes(iFEvol).varloc ==1))

        var_loc=tdata.FEvolumes(iFEvol).varloc;
        var_specified=1;
        
        %Wen Long, make sure when var_specified==1
        %data_packing is 0 (block) 
        %rule out conflict conditions (data_packing=1 (point) and var_loc=
        %1 (cell-center) cannot co-exist)
        %That is to say when var_loc=1, data_packing must be zero (block)
        %
        if(data_packing==1 && var_loc==1)
            s=-1;
            dispplay(['Error: datapacking of FEvolume ' int2str(iFEvol) ' (1 point) conflicts with varloc']);
            return
        end
   else
       warning(['var location of FEvolume ' int2str(iFEvol) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
       var_loc=0;
       var_specified=0;
       warning(['set default value as 0 (nodal)']);
   end
end

fwrite(fid_out,var_specified,'int32');


%
%check avaibility of x,y,z (coordinate variables) and size consistency
%

have_x=~isempty(find(strcmp(FEvolfields,'x')==1));
have_y=~isempty(find(strcmp(FEvolfields,'y')==1));
have_z=~isempty(find(strcmp(FEvolfields,'z')==1));
if(have_x && isempty(tdata.FEvolumes(iFEvol).x))
   have_x=0;
end
if(have_y && isempty(tdata.FEvolumes(iFEvol).y))
   have_y=0;
end
if(have_z && isempty(tdata.FEvolumes(iFEvol).y))
   have_z=0;
end

%volume is defined on 3D coordinates (x,y,z)
%v=v(x,y,z)
%(x,y,z) are coordinate vairables, must be nodal
%and exist
            
if(have_x && have_y && have_z)                
    %make sure x, y and z are vectors and have same size
    if(        isvector(tdata.FEvolumes(iFEvol).x)  ...
            && isvector(tdata.FEvolumes(iFEvol).y)  ...
            && isvector(tdata.FEvolumes(iFEvol).z)  ...
            && length(tdata.FEvolumes(iFEvol).x)== ...
               length(tdata.FEvolumes(iFEvol).y) ...
            && length(tdata.FEvolumes(iFEvol).x)== ...
               length(tdata.FEvolumes(iFEvol).z))
        NumNodes=length(tdata.FEvolumes(iFEvol).x);
    else
        s=-1;
        display(['Error: FEvolume ' int2str(iFEvol)  ...
                 ' x, y and z must be 1D vectors and have same size']);
        return;
    end
else
    s=-1;
    display(['Error: FEvolume ' int2str(iFEvol) ' must have ' ...
             ' x, y and z provided on all nodes ']);
    return;
end

 if(var_specified==1)
     

%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)

   %
   %give location for each variable in this zone
   %
   NV=tdata.Nvar;
   for iv=1:NV %NV is number of variables
       if(tdata.FEvolumes(iFEvol).varloc~=1)
           var_location = 0;  %nodal 
       else    %variable iv is cell centered
           var_location = 1; %cell centered
       end

       %make sure coordinate variables are nodal
       %assuming first 3 variables are coordinates for this volume!
       if(iv==1||iv==2||iv==3)  %x, y, z must be nodal
          var_location = 0; %Node
       end
       fwrite(fid_out,var_location,'int32');      
   end
 end

%
% are raw local 1-to-1 face neighbors supplied 
%
%         +-----------+
%         | INT32     |       Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+       These raw values are a compact form of the local 1-to-1 face
%                             neighbors are fully specified and therefore it will not
%                             perform auto face neighbor assignment, thereby improving
%                             Tecplot's time to first plot.
%                             See data section below for format details. ORDERED and
%                             FELINESEG zones must specify 0 for this value since raw
%                             face neighbors are not defined for these zone types.%

% Here we set it to false for ordered data 
% flase = 0 (int32), true = 1 (int32)
%
raw_1to1_face_supplied=0; %raw local 1-to-1 face neighbors not supplied
fwrite(fid_out,raw_1to1_face_supplied,'int32');

% For face connection (for each face in curren tzone), Tecplot takes the following
% *****The combination of cell and face numbers in the current zone must be unique*****
% ***** multiple entries are not allowed                                          *****
%-----------------------------------------------------------------------------------------------------
% FaceNeighbor Mode |  #Values  | Data
%-----------------------------------------------------------------------------------------------------
% LocalOneToOne     |   3       |  cz,fz,cz   (cz--cell number in current zone
%                   |           |              fz--cell face number in current zone
%                   |           |              cz--cell number of the neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% LocalOneToMany    |   nz+4    |  cz,fz,oz,nz,cz1,cz2,...,czn
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--number of cell face in current zone
%                   |           |        oz--face obscuration flag
%                   |           |        nz--number of neighboring cells for one-to-many options
%                   |           |        cz1--cell number of the neighbor cell in current zone
%                   |           |        cz2--cell number of 2nd neighbor cell in current zone
%                   |           |             ...
%                   |           |        czn--cell number of the n'th neighbor cell in current zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToOne    |   4       |  cz,fz,ZZ,CZ
%                   |           |       (cz--cell number in current zone
%                   |           |        fz--face number in current zone
%                   |           |        ZZ--remote zone number
%                   |           |        CZ--cell number of neighboring cell in remote zone)
%-----------------------------------------------------------------------------------------------------
% GlobalOneToMany   |   2*nz+4  |  cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%-----------------------------------------------------------------------------------------------------

%
% cz --cell in current zone
% fz --face of cell in current zone
% oz --face obscuration flag (only applies to one-to-many)
%         0 -- face partially obscured
%         1 -- face entirely obscured
% nz --number of cell or zone/cell associateions
% ZZ --remote Zone
% CZ --cell in remote zone
%
% cz,fz combinations must be unique. Additionally,Tecplot assume sthat with
% the one-to-one face neighbor modes a supplied cell face is entirely
% obscured by its neighbor. With one-to-many, the obscuration flag must be 
% supplied. 
%

% No. of miscellaneous user defined face neighbor connections
%         +-----------+
%         | INT32     |       Number of miscellaneous user defined face neighbor connections 
%         +-----------+       (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.

NoOfUserDefinedNeighbourConn = 0; % 0 for no, >0 for yes
fwrite(fid_out,NoOfUserDefinedNeighbourConn,'int32');

if (NoOfUserDefinedNeighbourConn ~=0)  
    
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)

     %Pls modify based on data
     user_defined_face_neighbor_mode= 0; %Local 1-to-1
     
     
     fwrite(fid_out,user_defined_face_neighbor_mode,'int32');
     
     if(zone_type>0)  %for non-oredered zones only
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
        
        %pls modify this value based on your data
        fe_face_complete_by_misc_face_neighbors=1;
        
        fwrite(fid_out,fe_face_complete_by_misc_face_neighbor,'int32');
     end
    
    
end

if(zone_type==0)  %for ordered zone

%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+

   %
   % ordered zone var 3 int32s for num points i -j -k
   %
   
%    Imax=
%    Jmax=
%    Kmax=
%    fwrite(fid_out,Imax,'int32');
%    fwrite(fid_out,Jmax,'int32');
%    fwrite(fid_out,Kmax,'int32');

else  %for non-ordered (FE) zones
   
%         if FE Zone:
%           +-----------+
%           | INT32     |     NumPts
%           +-----------+
% 		    if ZoneType is FEPOLYGON or FEPOLYHEDRON: 
% 		        +-----------+ 
% 			    | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					            no neighboring element. 
% 			    +-----------+ 
% 			    | INT32	    |   Total number of boundary connections. 
% 			    +-----------+
%
%           +-----------+
%           | INT32     |     NumElements.
%           +-----------+       
%           +-----------+
%           | INT32*3   |     ICellDim,JCellDim,KCellDim (for future use; set to zero)
%           +-----------+       

     % %uncomment and complete here

      NumPts=NumNodes  ; %give number of nodes
      fwrite(fid_out,NumPts,'int32');
      
         if(zone_type==6 || zone_type==7)  %FEPOLYGON or FEPOLYHEDRON
% 		        +-----------+ 
% 		        | INT32     |   NumFaces 
% 	    		+-----------+ 
% 			    +-----------+ 
% 			    | INT32     |   Total number of face nodes. For FEPOLYGON 
% 			    +-----------+   zones, this is NumFaces*2.  
%                                           For FEPOLYHEDRON, it is summation of all 
%                                           nodes that happen to be on a surface. If a
%                                           node belongs to 3 faces, it should be
%                                           counted 3 times.
% 			    +-----------+ 
% 			    | INT32	|   Total number of boundary faces. If any 
% 			    +-----------+   boundary faces exist, include one to represent 
% 					    no neighboring element. 
% 			    +-----------+ 
% 			    | INT32     |   Total number of boundary connections. 
% 			    +-----------+      

%          %uncomment the following and complete
%          %
%            %     
%            %Here we should use e2n structure's Faces,BoundaryFaces 
%            %etc members/properties array to to calculate the following:
%            %
%                 NumFaces=
%                 NumFaceNodes=
%                 NumBryFaces=
%                 NumBryConnections=
%                 fwrite(fid_out,NumFaces,'int32');
%                 fwrite(fid_out,NumFaceNodes,'int32');
%                 fwrite(fid_out,NumBryFaces,'int32');
%                 fwrite(fid_out,NumBryConnections,'int32');

         end
      
      %NumElements already calculated above. 
      
      fwrite(fid_out,NumElements,'int32');
      
      ICellDim=0;
      JCellDim=0;
      KCellDim=0;
      fwrite(fid_out,ICellDim,'int32');
      fwrite(fid_out,JCellDim,'int32');
      fwrite(fid_out,KCellDim,'int32');
    
end

% Zone auxiliary data pairs
%
%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%
% Auxiliary data may be used by text, macros, equations (if numeric) and
% add-ons in tecplot. It may be viewed directly in the Aux Data Page
% of the "Data Set Information" dialog (in "Data" menu)
% It must contain a name and value pair. The name must be a null-terminated
% character string and cannot contain spaces. the value must be a null
% terminated character string as well
%

%check if we have aux data for this FEvolume, true if both auxname and 
%auxval exist and are 'cell arrays or string' and neither of them is empty
have_auxdata=(   ~isempty(find(strcmp(FEvolfields,'auxname')==1))  ...
              && ~isempty(find(strcmp(FEvolfields,'auxval')==1)));   
if(have_auxdata)
   aux_data=1;  %with aux data. And we will retrieve each auxname and 
                %auxval pair here
    aux_names=tdata.FEvolumes(iFEvol).auxname;  %cell array or string
   aux_values=tdata.FEvolumes(iFEvol).auxval;   %cell array or string
   %
   %check if they are cell arrays (1D) or just a single string
   %
   if(iscellstr(aux_names) && iscellstr(aux_values))
      %or check if they are string (not array), just single string value
      N_aux=min(length(aux_names),length(aux_values));
                             %number of aux is obtained from aux_names
                             %and aux_values. 
      if(length(aux_names) > N_aux || length(aux_values)> N_aux)
          %give warnging if they do not match in size
          warning(['Warning: FEvolume ' int2str(iFEvol) ' auxname and auxval have inconsistent sizes']);
          warning(['Warning: only taking the first ' int2str(N_aux) ' values, ignoring the rest']);
      end
   else
       if(ischar(aux_names) && ischar(aux_values))
           N_aux=1;  %singe string given for aux_names and aux_values
           aux_names={aux_names};    %covnert string to cell
           aux_values={aux_values};  %convert string to cell
       else
           N_aux=0;  
       end
   end

   %if N_aux =0, then there is error, we can not have 
   %inconsistent auxname and auxval, they have to be both
   %cell array of strings or both a string value
   %give warning and quit
   
   if(N_aux==0)
       warning(['Warning: FEvolume ' int2str(iFEvol) ' auxname and auxval inconsistent']);
       warning(['set default as not to have aux data for this FEvolume']);
       aux_data=0;
   end
   
else
   aux_data=0;  %no aux data
end

%Wen Long: debug, do we have aux_data repeatedly provided
%for each aux pair when aux_data == 1? I think so.
if(aux_data==0)
    fwrite(fid_out,aux_data,'int32');  
end

if(aux_data==1)
    %repeat for each aux variable
    for iaux=1:N_aux
            
          fwrite(fid_out,aux_data,'int32');  %repeat for each aux variable
          
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      

          %
          %modify aux_name, aux_vale for each iaux
          %
          aux_name=aux_names{iaux};    
          aux_value_format=0;          %curently only allow string (0)
          aux_value=aux_values{iaux};  %currently only allow string values

          %make suer they are not empty, if yes, default to 'UndefinedAuxName1'
          %'UndefinedAuxName2',...
          %and 'UndefinedAuxVal1','UndefinedAuxVal2',...
          if(isempty(aux_name))
              aux_name=['UndefinedName' int2str(iaux)];
          end
          if(isempty(aux_value))
              aux_value=['UndefinedVal' int2str(iaux)];
          end
          
          %finally output this aux data pair
          plt_write_string(fid_out, aux_name);      %first string
          fwrite(fid_out,aux_value_format,'int32'); %second format
          plt_write_string(fid_out, aux_value);     %third value
          
    end
end
end

%--------

%
% Header section  v. Geometries
%

have_geoms=~isempty(find(strcmp(tdatanames,'geometries')==1));  %check if have geometries
%get number of geometries in tdata
if(have_geoms)  %make sure have geometries
    if(isstruct(tdata.geometries))  %make sure tdata.geometries is structure array
                                  %or at least structure 
       Ngeoms =length(tdata.geometries);
    else
       Ngeoms = 0;
    end
else
    Ngeoms=0;   
end

%loop though all geometries if any

for igeom=1:Ngeoms
    
    geomfields=fieldnames(tdata.geometries(igeom));  %find fields in
                                                     %this geometry

% special stuff such as circles, triangles, lines etc
% we do not have to have it for coordinate data, or mesh data
% looks like we can get away without it
%
%          +-----------+
%          | FLOAT32   |       Geometry marker.  Value = 399.0
%          +-----------+

    dummy_float32 = single(399.0);
    fwrite(fid_out,dummy_float32,'float32');

%          +-----------+
%          | INT32     |       Position CoordSys 0=Grid, 1=Frame, 2=FrameOffset(not used),
%          +-----------+                         3= OldWindow(not used), 4=Grid3D(New to V10)

    %check if there is 'cs' for this geometry, if yes output, if not
    %default to 0
    have_cs=~isempty(find(strcmp(geomfields,'cs')==1)); 
    if(have_cs && isempty(tdata.geometries(igeom).cs))
        have_cs=0;
    end
    if(~have_cs)
        geom_cs=0; 
        warning(['Geometry ' int2str(igeom) ' does not have cs (coordinate system)']);
        warning(['Set to default value as ' int2str(geom_cs)]);
    else
        geom_cs=tdata.geometries(igeom).cs;
    end
    
    fwrite(fid_out,geom_cs,'int32'); 

%          +-----------+
%          | INT32     |       Scope 0=Global 1=Local
%          +-----------+

    %check if there is 'scope' for this geometry, if yes output, if not
    %default to 0 (Global)
    have_scope=~isempty(find(strcmp(geomfields,'scope')==1)); 
    if(have_scope && isempty(tdata.geometries(igeom).scope))
        have_scope=0;
    end
    if(~have_scope)
        geom_scope=0;
        warning(['Geometry ' int2str(igeom) ' does not have scope (0-global, 1-local)']);
        warning(['Set to default value as ' int2str(geom_scope)]);
    else
        geom_scope=tdata.geometries(igeom).scope;
    end
    
    fwrite(fid_out,geom_scope,'int32'); 

%          +-----------+
%          | INT32     |       DrawOrder 0=After, 1=Before
%          +-----------+

    %check if there is 'draworder' for this geometry, if yes output, if not
    %default to 0 (After)
    have_draworder=~isempty(find(strcmp(geomfields,'draworder')==1)); 
    if(have_draworder && isempty(tdata.geometries(igeom).draworder))
        have_draworder=0;
    end
    if(~have_draworder)
        geom_draworder=0; 
        warning(['Geometry ' int2str(igeom) ' does not have draworder (0-global, 1-local)']);
        warning(['Set to default value as ' int2str(geom_draworder)]);
    else
        geom_draworder=tdata.geometries(igeom).draworder;
    end
    fwrite(fid_out,geom_draworder,'int32'); 

%          +-----------+
%          | FLOAT64*3 |       (X or Theta),(Y or R),(Z or dummy)  
%          +-----------+       i.e. the starting location
     %check if have x or theta, y or r, z 
     %(must have x or theta, y or r.
     % If no z, z will default to zero
     % If both x and theta are available , take x
     % If both y and r are available, take y
     have_x=~isempty(find(strcmp(geomfields,'x')==1));
     have_theta=~isempty(find(strcmp(geomfields,'theta')==1));
     have_y=~isempty(find(strcmp(geomfields,'y')==1));
     have_r=~isempty(find(strcmp(geomfields,'r')==1));
     have_z=~isempty(find(strcmp(geomfields,'z')==1));

     if(have_x && isempty(tdata.geometries(igeom).x))
         have_x=0;
     end
     if(have_y && isempty(tdata.geometries(igeom).y))
         have_y=0;
     end
     if(have_z && isempty(tdata.geometries(igeom).z))
         have_z=0;
     end
     if(have_theta && isempty(tdata.geometries(igeom).theta))
         have_theta=0;
     end
     if(have_r && isempty(tdata.geometries(igeom).r))
         have_r=0;
     end

     if(~have_x && ~have_theta)
         s=-1;
         display(['Error: geometry ' int2str(igeom) ' does not have ' ...
                  'x or theta']);
         return           
     else
         if(have_theta)
             geom_x=tdata.geometries(igeom).theta;
         end
         if(have_x)   %this is after theta, so x will over ride theta when provided
             geom_x=tdata.geometries(igeom).x;
         end
     end
     if(~have_y && ~have_r)
         s=-1;
         display(['Error: geometry ' int2str(igeom) ' does not have ' ...
                  'y or r']);
         return                    
     else
         if(have_r)
             geom_y=tdata.geometries(igeom).r;
         end
         if(have_y)  %this is after r, so y will over rider r when provided
             geom_y=tdata.geometries(igeom).y;
         end
     end
     if(~have_z)
        geom_z=0;  %default to zero
        warning(['Geometry ' int2str(igeom) ' does not have z ']);
        warning(['Set to default as ' num2str(geom_z)]);
     else
        geom_z=tdata.geometries(igeom).z;
     end
     %make sure x and y exist at same time, or do not exist at same time
     if(xor(have_x,have_y))
         s=-1;
         display(['Error: geometry ' int2str(igeom) ' must have x and y ' ...
                  '       at the same time or neither of them']);
         return
     end
     %make sure theta and r exist at same time, or do not exist at same
     %time
     if(xor(have_r,have_theta))
         s=-1;
         display(['Error: geometry ' int2str(igeom) ' must have r and theta ' ...
                  '       at the same time or neither of them']);
         return
     end

     %make sure anchoring location coordinate is from 0 to 1 
     %when coordinate system is Frame (1)
     if(geom_cs==1)
        if(geom_x >=1 || geom_x <0)
          display(['Warning: geometry ' int2str(igeom) ' x or r anchoring should be between 0 and 1 when cs is 1 (Frame)']);
        end
        if(geom_y >=1 || geom_y <0)
          display(['Warning: geometry ' int2str(igeom) ' y or theta anchoring should be between 0 and 1 when cs is 1 (Frame)']);
        end
     end

     fwrite(fid_out,geom_x,'float64'); 
     fwrite(fid_out,geom_y,'float64'); 
     fwrite(fid_out,geom_z,'float64'); 
     
%          +-----------+
%          | INT32     |       Zone (0=all)
%          +-----------+

    %check if there is 'zone' for this geometry, if yes output, if not
    %default to 0 (all)
    have_zone=~isempty(find(strcmp(geomfields,'zone')==1)); 
    if(have_zone && isempty(tdata.geometries(igeom).zone))
        have_zone=0;
    end
    if(~have_zone)
        geom_zone=0; 
        warning(['Geometry ' int2str(igeom) ' does not have specified zone number']);
        warning(['Set to default value as ' int2str(geom_zone)]);
    else
        geom_zone=tdata.geometries(igeom).zone;
    end

    fwrite(fid_out,geom_zone,'int32'); 
    
%          +-----------+
%          | INT32     |       Color
%          +-----------+

    %check if there is 'color' for this geometry, if yes output, if not
    %default to 0
    have_color=~isempty(find(strcmp(geomfields,'color')==1)); 
    if(have_color && isempty(tdata.geometries(igeom).color))
        have_color=0;
    end
    if(~have_color)
        geom_color=0; 
        warning(['Geometry ' int2str(igeom) ' does not have specified color']);
        warning(['Set to default value as ' int2str(geom_color)]);
    else
        geom_color=tdata.geometries(igeom).color;
    end

    fwrite(fid_out,geom_color,'int32'); 

%          +-----------+
%          | INT32     |       FillColor
%          +-----------+


    %check if there is 'fillcolor' for this geometry, if yes output, if not
    %default to 0 
    have_fillcolor=~isempty(find(strcmp(geomfields,'fillcolor')==1)); 
    if(have_fillcolor && isempty(tdata.geometries(igeom).fillcolor))
        have_fillcolor=0;
    end
    if(~have_fillcolor)
        geom_fillcolor=0; 
        warning(['Geometry ' int2str(igeom) ' does not have specified fillcolor']);
        warning(['Set to default value as ' int2str(geom_fillcolor)]);
    else
        geom_fillcolor=tdata.geometries(igeom).fillcolor;
    end

    fwrite(fid_out,geom_fillcolor,'int32'); 

%          +-----------+
%          | INT32     |       IsFilled (0=no 1=yes)
%          +-----------+

    %check if there is 'isfilled' for this geometry, if yes output, if not
    %default to 0 (no)
    have_isfilled=~isempty(find(strcmp(geomfields,'isfilled')==1)); 
    if(have_isfilled && isempty(tdata.geometries(igeom).isfilled))
        have_isfilled=0;
    end
    if(~have_isfilled)
        geom_isfilled=0; 
        warning(['Geometry ' int2str(igeom) ' does not have isfilled flag']);
        warning(['Set to default value as ' int2str(geom_isfilled)]);
    else
        geom_isfilled=tdata.geometries(igeom).isfilled;
    end
    fwrite(fid_out,geom_isfilled,'int32'); 

%          +-----------+
%          | INT32     |       GeomType  0=Line,  1=Rectangle, 2=Square,
%          +-----------+                 3=Circle,4=ellipse,   5=Line3D
%

    have_geomtype=~isempty(find(strcmp(geomfields,'geomtype')==1)); 
    if(have_geomtype && isempty(tdata.geometries(igeom).geomtype))
        have_geomtype=0;
    end
    if(~have_geomtype)
        s=-1;
        display(['Error: Geometry ' int2str(igeom) ' does not have geomtype']);
        display(['       Set to default value as ' int2str(geom_isfilled)]);
        return;
    else
        geom_type=tdata.geometries(igeom).geomtype;
    end
    if(geom_type<0 || geom_type>5)
        s=-1;
        display(['Error: Geometry ' int2str(igeom) ' geomtype must be 0,1,2,3,4 or 5']);
        return;
    end


    fwrite(fid_out,geom_type,'int32'); 


%          +-----------+
%          | INT32     |       LinePattern  0=Solid 1=Dashed 2=DashDot 
%          +-----------+                    3=DashDotDot,4=Dotted
%                                           5=LongDash

    have_linepattern=~isempty(find(strcmp(geomfields,'linepattern')==1));
    if(have_linepattern && isempty(tdata.geometries(igeom).linepattern))
        have_linepattern=0;
    end
    if(~have_linepattern)
        geom_linepattern=0;
        warning(['Geometry ' int2str(igeom) ' does not have linepattern']);
        warning(['Set to default value as ' int2str(geom_linepattern)]);
    else
        geom_linepattern=tdata.geometries(igeom).linepattern;
    end

    fwrite(fid_out,geom_linepattern,'int32');

%          +-----------+ 
%          | FLOAT64   |       Pattern Length (in frame unit)
%          +-----------+
    
    have_patternlength=~isempty(find(strcmp(geomfields,'patternlength')==1));
    if(have_patternlength && isempty(tdata.geometries(igeom).patternlength))
        have_patternlength=0;
    end
    if(~have_patternlength)
        geom_patternlength=0.01;
        warning(['Geometry ' int2str(igeom) ' does not have patternlength']);
        warning(['Set to default value as ' int2str(geom_patternlength)]);
    else
        geom_patternlength=tdata.geometries(igeom).patternlength;
    end

    fwrite(fid_out,geom_patternlength,'float64');


%          +-----------+
%          | FLOAT64   |       Line Thickness
%          +-----------+
    have_linethickness=~isempty(find(strcmp(geomfields,'linethickness')==1));
    if(have_linethickness && isempty(tdata.geometries(igeom).linethickness))
        have_linethickness=0;
    end
    if(~have_linethickness)
        geom_linethickness=0.001;  %this is 0.1%
        warning(['Geometry ' int2str(igeom) ' does not have linethickness']);
        warning(['Set to default value as ' int2str(geom_linethickness)]);
    else
        geom_linethickness=tdata.geometries(igeom).linethickness;
    end
    fwrite(fid_out,geom_linethickness,'float64');


%          +-----------+
%          | INT32     |  NumEllipsePts
%          +-----------+   Number of points used to approximate
%                          Ellipses and circles. Default is 72

    have_numellipsepts=~isempty(find(strcmp(geomfields,'numellipsepts')==1));
    if(have_numellipsepts && isempty(tdata.geometries(igeom).numellipsepts))
        have_numellipsepts=0;
    end
    if(~have_numellipsepts)
        geom_numellipsepts=72;
        warning(['Geometry ' int2str(igeom) ' does not have numellipsepts']);
        warning(['Set to default value as ' int2str(geom_numellipsepts)]);
    else
        geom_numellipsepts=tdata.geometries(igeom).numellipsepts;
    end

    fwrite(fid_out,geom_numellipsepts,'int32');

%                            geometries.arrowheadstyle  ---0-plain
%                                                          1-hollow
%                                                          2-filled
%                            geometries.arrowheadattach ---0-none
%                                                          1-beginning
%                                                          2-end
%                                                          3-both
%                            geometries.arrowheadsize   ---size of arrow
%                                                          head in frame units
%                            geometries.arrowheadangle  ---angle of arrow head
%                                                          in degrees


%          +-----------+
%          | INT32     |  Arrowhead Style 0=Plain, 1=Filled, 2=Hollow
%          +-----------+
    have_arrowheadstyle=~isempty(find(strcmp(geomfields,'arrowheadstyle')==1));
    if(have_arrowheadstyle && isempty(tdata.geometries(igeom).arrowheadstyle))
        have_arrowheadstyle=0;
    end

    if(~have_arrowheadstyle)
        geom_arrowheadstyle=0;
        warning(['Geometry ' int2str(igeom) ' does not have arrowheadstyle']);
        warning(['Set to default value as ' int2str(geom_arrowheadstyle)]);
    else
        geom_arrowheadstyle=tdata.geometries(igeom).arrowheadstyle;
    end

    fwrite(fid_out,geom_arrowheadstyle,'int32');



%          +-----------+
%          | INT32     |  Arrowhead Attachment 0=None, 1=Beg, 2=End, 3=Both
%          +-----------+

    have_arrowheadattach=~isempty(find(strcmp(geomfields,'arrowheadattach')==1));
    if(have_arrowheadattach && isempty(tdata.geometries(igeom).arrowheadattach))
        have_arrowheadattach=0;
    end

    if(~have_arrowheadattach)
        geom_arrowheadattach=0;
        warning(['Geometry ' int2str(igeom) ' does not have arrowheadattach']);
        warning(['Set to default value as ' int2str(geom_arrowheadattach)]);
    else
        geom_arrowheadattach=tdata.geometries(igeom).arrowheadattach;
    end

    fwrite(fid_out,geom_arrowheadattach,'int32');

%          +-----------+
%          | FLOAT64   |  Arrowhead Size (in frame units)
%          +-----------+

    have_arrowheadsize=~isempty(find(strcmp(geomfields,'arrowheadsize')==1));
    if(have_arrowheadsize && isempty(tdata.geometries(igeom).arrowheadsize))
        have_arrowheadsize=0;
    end
    if(~have_arrowheadsize)
        geom_arrowheadsize=0.01;
        warning(['Geometry ' int2str(igeom) ' does not have arrowheadsize']);
        warning(['Set to default value as ' int2str(geom_arrowheadsize)]);
    else
        geom_arrowheadsize=tdata.geometries(igeom).arrowheadsize;
    end

    fwrite(fid_out,geom_arrowheadsize,'int32');

%          +-----------+
%          | FLOAT64   |  Arrowhead Angle
%          +-----------+

    have_arrowheadangle=~isempty(find(strcmp(geomfields,'arrowheadangle')==1));
    if(have_arrowheadangle && isempty(tdata.geometries(igeom).arrowheadangle))
        have_arrowheadangle=0;
    end
    if(~have_arrowheadangle)
        geom_arrowheadangle=0;
        warning(['Geometry ' int2str(igeom) ' does not have arrowheadangle']);
        warning(['Set to default value as ' int2str(geom_arrowheadangle)]);
    else
        geom_arrowheadangle=tdata.geometries(igeom).arrowheadangle;
    end

    fwrite(fid_out,geom_arrowheadangle,'int32');


%          +-----------+
%          | IN32*N    |  Macro Function Command (string: N = Length+1)
%          +-----------+

    %this is not supported here. Simply give no commands
    geom_MFC_dummy='none';
    plt_write_string(fid_out, geom_MFC_dummy);


%          +-----------+
%          | INT32     |  Polyline Field Data Type 1=Float, 2=Double  (GTYPE)
%          +-----------+

    have_datatype=~isempty(find(strcmp(geomfields,'datatype')==1));
    if(have_datatype && isempty(tdata.geometries(igeom).datatype))
        have_datatype=0;
    end
    if(~have_datatype)
        geom_datatype=1;
        warning(['Geometry ' int2str(igeom) ' does not have datatype']);
        warning(['Set to default value as ' int2str(geom_datatype)]);
    else
        geom_datatype=tdata.geometries(igeom).datatype;
    end

    if(geom_datatype ~=1 && geom_datatype ~=2)
       s=-1;
       display(['Error, geometry ' int2str(igeom) ' datatype must be 1(Float) or 2(Double)']);
       return;
    end

    fwrite(fid_out,geom_datatype,'int32');


%          +-----------+
%          | INT32     |  Clipping (0=ClipToAxes,1=ClipToViewport,2=ClipToFrame)
%          +-----------+
% 

    have_clipping=~isempty(find(strcmp(geomfields,'clipping')==1));
    if(have_clipping && isempty(tdata.geometries(igeom).clipping))
        have_clipping=0; 
    end
    if(~have_clipping)
        geom_clipping=0;
        warning(['Geometry ' int2str(igeom) ' does not have clipping']);
        warning(['Set to default value as ' int2str(geom_clipping)]);
    else
        geom_clipping=tdata.geometries(igeom).clipping;
    end
    fwrite(fid_out,geom_clipping,'int32');

    %provide data describing the geometry based on geom_type
%
%    geometries.data  ---data array that describes
%                        the geometry
%                        for SQUARE    - one number of
%                                        side length
%                        for RECTANGLE - two numbers of
%                                        width and height
%                        for CIRCLE    - one number of radius
%                        for ELLIPSE   - two numbers of hor-
%                                        izontal axis and vertical
%                                        axis length
%
%                        for LINE or LINE3D, the data is a
%                                    colloection of up to 50 polylines
%                                    and each polyline is defined by
%                                    a number of XY or XYZ coordinates
%                                    The cooridnates are relative
%                                    to anchor point.

   %makke sure data is available

   have_data=~isempty(find(strcmp(geomfields,'data')==1));
   if(have_data && isempty(tdata.geometries(igeom).data))
        have_data=0 ;
   end
   if(~have_data)
      s=-1;
      display(['Error: geometry ' int2str(igeom) ' does not have data !']);
      return;
   end

   if(geom_type==0 ||geom_type==5)  %Line or Line3D
%         GeomType can be  0=Line, 1=Rectangle 2=Square,
%                          3=Circle, 4=ellipse,5=Line3D
%

             %make sure data is a structure storing 
             %polylines 
             %
             %    data(iline).x -- x coordinate vector
             %    data(iline).y -- y coordinate vector
             %    data(ilien).z -- z coordinate vector
             %

           num_polylines=0;
           is_data_struct=isstruct(tdata.geometries(igeom).data);
           if(~is_data_struct)
              s=-1;
              display(['Error: geometry ' int2str(igeom) ' data should be a structure vector']);
              display(['       storing (x,y) for Line or (x,y,z) for Line3D']);
              return;
           else  
             num_polylines=length(tdata.geometries(igeom).data);
             
           end
           %make sure there is at least one polyline
           if(num_polylines<1)
              s=-1;
              display(['Error: geometry ' int2str(igeom) ' data must contain at least one set']);
              display(['       of polyline']);
              return
           end

%            If the geometry type is line or line 3D then:
%                +-----------+
%                | INT32     |       Number of polylines
%                +-----------+

               fwrite(fid_out,num_polylines,'int32');

               for iline=1:num_polylines          %Loop through all polylines
                   polylinefields=fieldnames(tdata.geometries(igeom).data(iline));
                  
                   have_poly_x=~isempty(find(strcmp(polylinefields,'x')==1)); 
                   have_poly_y=~isempty(find(strcmp(polylinefields,'y')==1));
                   have_poly_z=~isempty(find(strcmp(polylinefields,'z')==1));
                   if(have_poly_x && isempty(tdata.geometries(igeom).data(iline).x))
                      have_poly_x=0;
                   end
                   if(have_poly_y && isempty(tdata.geometries(igeom).data(iline).y))
                      have_poly_y=0;
                   end
                   if(have_poly_z && isempty(tdata.geometries(igeom).data(iline).z))
                      have_poly_z=0;
                   end


                   if(~have_poly_x||~have_poly_y)
                       s=-1;
                       display(['Error: geometry ' int2str(igeom) ' polyline ' ...
                                        int2str(iline) ' must have both x and y']);
                       return
                   end

                   %make sure x and y are vectors (row or column)
                   is_poly_x_vector=isvector(tdata.geometries(igeom).data(iline).x);
                   is_poly_y_vector=isvector(tdata.geometries(igeom).data(iline).y);

                   NumPts_x=0;
                   NumPts_y=0;
                   if(is_poly_x_vector && is_poly_y_vector)
                      NumPts_x=length(tdata.geometries(igeom).data(iline).x);
                      NumPts_y=length(tdata.geometries(igeom).data(iline).y);
                   else
                      s=-1;
                      display(['Error: geometry ' int2str(igeom) ' polyline ' ...
                                       int2str(iline) ' x and y must be a vector']);
                      return;
                   end

                   %make sure there are at least two points to define a line
                   if(NumPts_x<2 || NumPts_y<2 || NumPts_x~=NumPts_y)
                     s=-1;
                     display(['Error: geometry ' int2str(igeom) ' polyline ']);
                     display(['       ' int2str(iline) ' x and y have at least two points']);
                     display(['         and same number of points']);
                   end

%
%                +-----------+
%                | INT32     |       Number of points, line 1.
%                +-----------+

                   %output Number of points for this polyline
                   fwrite(fid_out,NumPts_x,'int32');

%                +-----------+
%                | GTYPE*N   |       X-block geometry points N=NumPts
%                +-----------+
                   for iPoint=1:NumPts_x
                      switch(geom_datatype)
                        case 1  %float
                           fwrite(fid_out,tdata.geometries(igeom).data(iline).x(iPoint),'float32');
                        case 2  %double
                           fwrite(fid_out,tdata.geometries(igeom).data(iline).x(iPoint),'float64');
                      end
                   end

%                +-----------+
%                | GTYPE*N   |       Y-block geometry points N=NumPts
%                +-----------+
                   for iPoint=1:NumPts_y
                      switch(geom_datatype)
                        case 1  %float
                           fwrite(fid_out,tdata.geometries(igeom).data(iline).y(iPoint),'float32');
                        case 2  %double
                           fwrite(fid_out,tdata.geometries(igeom).data(iline).y(iPoint),'float64');
                      end
                   end
 
                   if(geom_type==5) %Line3D (Grid3D Only)
                    %makke sure z is available
                    %also make sure coordinate system is Grid3D
                    if(~have_poly_z || geom_cs~=4)
                        s=-1;
                        display(['Error: geometry ' int2str(igeom) ' polyline ']);
                        display(['       ' int2str(iline) ' must have z vector ']);
                        display(['       and coordinate system (cs) must be 4(Grid3D)']);
                        display(['       to handle 3D polylines']);
                        return;
                    end

                    is_poly_z_vector=isvector(tdata.geometries(igeom).data(iline).z);
                    NumPts_z=0;
                    if(~is_poly_z_vector)
                        s=-1;
                        display(['Error: geometry ' int2str(igeom) ' polyline ' int2str(iline)]);
                        display(['       must have z as a vector (row or column)']);
                        return;
                    else
                        NumPts_z=length(tdata.geometries(igeom).data(iline).z);
                    end
                    if(NumPts_z~=NumPts_x||NumPts_z~=NumPts_y)
                        s=-1;
                        display(['Error: geometry ' int2str(igeom) ' polyline ' int2str(iline)]);
                        display(['       z must have same size as x and y']);
                        return
                    end

%                   +-----------+
%                   | GTYPE*N   |       Z-block geometry points N=NumPts (Grid3D Only)
%                   +-----------+

                      for iPoint=1:NumPts_z
                          switch(geom_datatype)
                           case 1  %float
                              fwrite(fid_out,tdata.geometries(igeom).data(iline).z(iPoint),'float32');
                           case 2  %double
                              fwrite(fid_out,tdata.geometries(igeom).data(iline).z(iPoint),'float64');
                          end
                       end
                   end
               end
   end

   if(geom_type==1)
%            If the geometry type is Rectangle then
%                +-----------+
%                | GTYPE*2   |       X and Y offset for far corner of rectangle
%                +-----------+

         is_data_vector=isvector(tdata.geometries(igeom).data);
         if(~is_data_vector)
           s=-1;
           display(['Error: geometry ' int2str(igeom) ' rectangle ']);
           display(['       data must be a vector giving X and Y offset']);
           return
         else
           geom_rect_x=tdata.geometries(igeom).data(1);
           geom_rect_y=tdata.geometries(igeom).data(2);
         end
         if(isfinite(geom_rect_x) && isfinite(geom_rect_y))
              switch(geom_datatype)
                 case 1  %float
                    fwrite(fid_out,geom_rect_x,'float32');
                    fwrite(fid_out,geom_rect_y,'float32');
                 case 2  %double
                    fwrite(fid_out,geom_rect_x,'float64');
                    fwrite(fid_out,geom_rect_y,'float64');
              end
         else
            s=-1;
            display(['Error: geometry ' int2str(igeom) ' rectangle' ]);
            display(['       data must have finite X and Y offset values']);
            return;
         end

   end

   if(geom_type==2)
%            If the geometry type is Square then
%                +-----------+
%                | GTYPE     |       Width
%                +-----------+
%

        is_data_numeric=isnumeric(tdata.geometries(igeom).data);
        if(~is_data_numeric)
           s=-1;
           display(['Error: geometry ' int2str(igeom) ' square ']);
           display(['       data must be a numeric value giving Width']);
           return;
        else
           geom_square_width=tdata.geometries(igeom).data(1);
        end

        if(isfinite(geom_square_width))
           switch(geom_datatype)
                 case 1  %float
                    fwrite(fid_out,geom_square_width,'float32');
                 case 2  %double
                    fwrite(fid_out,geom_square_width,'float64');
           end
        else
            s=-1;
            display(['Error: geometry ' int2str(igeom) ' square' ]);
            display(['       data must have finite value for Width']);
            return;
        end
   end

   if(geom_type==3)

%            If the geometry type is Circle then
%                +-----------+
%                | GTYPE     |       Radius
%                +-----------+
% 

        is_data_numeric=isnumeric(tdata.geometries(igeom).data);
        if(~is_data_numeric)
           s=-1;
           display(['Error: geometry ' int2str(igeom) ' cicle ']);
           display(['       data must be a numeric value giving Radius']);
           return;
        else
           geom_circle_radius=tdata.geometries(igeom).data(1);
        end
        if(isfinite(geom_circle_radius))
           switch(geom_datatype)
                 case 1  %float
                    fwrite(fid_out,geom_circle_radius,'float32');
                 case 2  %double
                    fwrite(fid_out,geom_circle_radius,'float64');
           end
        else
            s=-1;
            display(['Error: geometry ' int2str(igeom) ' circle' ]);
            display(['       data must have finite value for Radius']);
            return;
        end
   end

   if(geom_type==4)
%            If the geometry type is Ellipse then
%                +-----------+
%                | GTYPE*2   |       X and Y Radii (semi x and semi y axis)
%                +-----------+
% 

         is_data_vector=isvector(tdata.geometries(igeom).data);
         if(~is_data_vector)
           s=-1;
           display(['Error: geometry ' int2str(igeom) ' ellipse ']);
           display(['       data must be a vector giving semi-X and semi-Y axis length']);
           return
         else
           geom_ellipse_x=tdata.geometries(igeom).data(1);
           geom_ellipse_y=tdata.geometries(igeom).data(2);
         end
         if(isfinite(geom_ellipse_x) && isfinite(geom_ellipse_y))
              switch(geom_datatype)
                 case 1  %float
                    fwrite(fid_out,geom_ellipse_x,'float32');
                    fwrite(fid_out,geom_ellipse_y,'float32');
                 case 2  %double
                    fwrite(fid_out,geom_ellipse_x,'float64');
                    fwrite(fid_out,geom_ellipse_y,'float64');
              end
         else
            s=-1;
            display(['Error: geometry ' int2str(igeom) ' ellipse' ]);
            display(['       data must have finite semi-X and semi-Y axis length values']);
            return;
         end
   end

end  %end of igeom loop

%
% Header section vi. Text
%

have_texts=~isempty(find(strcmp(tdatanames,'texts')==1));  %check if have texts
%get number of texts in tdata
if(have_texts)  %make sure have texts
    if(isstruct(tdata.texts))  %make sure tdata.texts is structure array
                               %or at least structure 
       NTexts =length(tdata.texts);
    else
       NTexts = 0;
    end
else
    NTexts=0;   
end

%loop though all texts if any

for itxt=1:NTexts


    textfields=fieldnames(tdata.texts(itxt));  %find fields in this text


%          +-----------+
%          | FLOAT32   |       Text marker.  Value=499.0
%          +-----------+
     dummy_float32 = single(499.0);
     fwrite(fid_out,dummy_float32,'float32');


%          +-----------+
%          | INT32     |       Position CoordSys 0=Grid, 1=Frame, 
%          +-----------+                         2=FrameOffset(not used),
%                                                3= OldWindow(not used),
%                                                4=Grid3D(New to V10)

    %check if there is 'cs' for this text, if yes output, if not
    %default to 0
    have_cs=~isempty(find(strcmp(textfields,'cs')==1));
    if(have_cs && isempty(tdata.texts(itxt).cs))
        have_cs=0;
    end
    if(~have_cs)
        text_cs=0;
        warning(['text ' int2str(itxt) ' does not have cs (coordinate system)']);
        warning(['Set to default value as ' int2str(text_cs)]);
    else
        text_cs=tdata.texts(itxt).cs;
    end
    fwrite(fid_out,text_cs,'int32');


%          +-----------+     
%          | INT32     |       Scope 0=Global 1=Local
%          +-----------+

    %check if there is 'scope' for this text, if yes output, if not
    %default to 0 (Global)
    have_scope=~isempty(find(strcmp(textfields,'scope')==1));
    if(have_scope && isempty(tdata.texts(itxt).scope))
        have_scope=0;
    end
    if(~have_scope)
        text_scope=0;
        warning(['text ' int2str(itxt) ' does not have scope (0-global, 1-local)']);
        warning(['Set to default value as ' int2str(text_scope)]);
    else
        text_scope=tdata.texts(itxt).scope;
    end

    fwrite(fid_out,text_scope,'int32');


%          +-----------+
%          | FLOAT64*3 |       (X or Theta),(Y or R),(Z or dummy)
%          +-----------+        Starting Location

     %check if have x or theta, y or r, z
     %(must have x or theta, y or r.
     % If no z, z will default to zero
     % If both x and theta are available , take x
     % If both y and r are available, take y
     have_x=~isempty(find(strcmp(textfields,'x')==1));
     have_theta=~isempty(find(strcmp(textfields,'theta')==1));
     have_y=~isempty(find(strcmp(textfields,'y')==1));
     have_r=~isempty(find(strcmp(textfields,'r')==1));
     have_z=~isempty(find(strcmp(textfields,'z')==1));
     if(have_x && isempty(tdata.texts(itxt).x))
        have_x=0;
     end
     if(have_y && isempty(tdata.texts(itxt).y))
        have_y=0;
     end
     if(have_z && isempty(tdata.texts(itxt).z))
        have_z=0;
     end
     if(have_theta && isempty(tdata.texts(itxt).theta))
        have_theta=0;
     end
     if(have_r && isempty(tdata.texts(itxt).r))
        have_r=0;
     end

     if(~have_x && ~have_theta)
         s=-1;
         display(['Error: text ' int2str(itxt) ' does not have ' ...
                  'x or theta']);
         return
     else
         if(have_theta)
             text_x=tdata.texts(itxt).theta;
         end
         if(have_x)
             text_x=tdata.texts(itxt).x;
         end
     end
     if(~have_y && ~have_r)
         s=-1;
         display(['Error: text ' int2str(itxt) ' does not have ' ...
                  'y or r']);
         return
     else
         if(have_r)
             text_y=tdata.texts(itxt).r;
         end
         if(have_y)
             text_y=tdata.texts(itxt).y;
         end
     end
     if(~have_z)
        text_z=0;  %default to zero
        warning(['text ' int2str(itxt) ' does not have z ']);
        warning(['Set to default as ' num2str(text_z)]);
     else
        text_z=tdata.texts(itxt).z;
     end
     %make sure x and y exist at same time, or do not exist at same time
     if(xor(have_x,have_y))
         s=-1;
         display(['Error: text ' int2str(itxt) ' must have x and y ' ...
                  '       at the same time or neither of them']);
         return
     end
     %make sure theta and r exist at same time, or do not exist at same
     %time
     if(xor(have_r,have_theta))
         s=-1;
         display(['Error: text ' int2str(itxt) ' must have r and theta ' ...
                  '       at the same time or neither of them']);
         return
     end
     fwrite(fid_out,text_x,'float64');
     fwrite(fid_out,text_y,'float64');
     fwrite(fid_out,text_z,'float64');

%          +-----------+
%          | INT32     |       FontType
%          +-----------+

    %check if there is 'fonttype' for this text, if yes output, if not
    %default to 5 (TIMES)
    have_fonttype=~isempty(find(strcmp(textfields,'fonttype')==1));
    if(have_fonttype && isempty(tdata.texts(itxt).fonttype))
        have_fonttype=0;
    end
    if(~have_fonttype)
        text_fonttype=5; 
        warning(['text ' int2str(itxt) ' does not have fonttype']);
        warning(['Set to default value as ' int2str(text_fonttype)]);
    else
        text_fonttype=tdata.texts(itxt).fonttype;
    end

    fwrite(fid_out,text_fonttype,'int32');


%          +-----------+
%          | INT32     |       Character Height Units 0=Grid, 
%          +-----------+                              1=Frame, 2=Point

    %check if there is 'charheightunit' for this text, if yes output, if not
    %default to 0
    have_charheightunit=~isempty(find(strcmp(textfields,'charheightunit')==1));
    if(have_charheightunit && isempty(tdata.texts(itxt).charheightunit))
        have_charheightunit=0;
    end
    if(~have_charheightunit)
        text_charheightunit=2; %default is point
        warning(['text ' int2str(itxt) ' does not have charheightunit(0-Grid,1=Frame,2-Point)']);
        warning(['Set to default value as ' int2str(text_charheightunit)]);
    else
        text_charheightunit=tdata.texts(itxt).charheightunit;
    end
    fwrite(fid_out,text_charheightunit,'int32');


%          +-----------+
%          | FLOAT64   |       Height of characters
%          +-----------+

    %check if there is 'charheight' for this text, if yes output, if not
    %default to 12
    have_charheight=~isempty(find(strcmp(textfields,'charheight')==1));
    if(have_charheight && isempty(tdata.texts(itxt).charheight))
        have_charheight=0;
    end
    if(~have_charheight)
        text_charheight=12; %default is 12 point
        warning(['text ' int2str(itxt) ' does not have charheight']);
        warning(['Set to default value as ' int2str(text_charheight)]);
    else
        text_charheight=tdata.texts(itxt).charheight;
    end
    fwrite(fid_out,text_charheight,'float64');

%          +-----------+
%          | INT32     |       Text Box type 0=NoBox 1=Hollow 2=Filled
%          +-----------+

    %check if there is 'boxtype' for this text, if yes output, if not
    %default to 0
    have_boxtype=~isempty(find(strcmp(textfields,'boxtype')==1));
    if(have_boxtype && isempty(tdata.texts(itxt).boxtype))
        have_boxtype=0;
    end
    if(~have_boxtype)
        text_boxtype=0; %default is 12 point
        warning(['text ' int2str(itxt) ' does not have boxtype (0-NoBox,1-Hollow,2-Filled)']);
        warning(['Set to default value as ' int2str(text_boxtype)]);
    else
        text_boxtype=tdata.texts(itxt).boxtype;
    end
    fwrite(fid_out,text_boxtype,'int32');

%          +-----------+
%          | FLOAT64   |       Text Box Margin
%          +-----------+
    
    have_boxmargin=~isempty(find(strcmp(textfields,'boxmargin')==1));
    if(have_boxmargin && isempty(tdata.texts(itxt).boxmargin))
        have_boxmargin=0;
    end
    if(~have_boxmargin)
        text_boxmargin=1; %default is 1
        warning(['text ' int2str(itxt) ' does not have boxmargin']);
        warning(['Set to default value as ' int2str(text_boxmargin)]);
    else
        text_boxmargin=tdata.texts(itxt).boxmargin;
    end
    fwrite(fid_out,text_boxmargin,'float64');


%          +-----------+
%          | FLOAT64   |       Text Box Margin Linewidth
%          +-----------+

    have_boxlinewidth=~isempty(find(strcmp(textfields,'boxlinewidth')==1));
    if(have_boxlinewidth && isempty(tdata.texts(itxt).boxlinewidth))
        have_boxlinewidth=0;
    end
    if(~have_boxlinewidth)
        text_boxlinewidth=1;  %default is 1
        warning(['text ' int2str(itxt) ' does not have boxlinewidth']);
        warning(['Set to default value as ' int2str(text_boxlinewidth)]);
    else
        text_boxlinewidth=tdata.texts(itxt).boxlinewidth;
    end
    fwrite(fid_out,text_boxlinewidth,'float64');

%          +-----------+
%          | INT32     |       Text Box Outline Color
%          +-----------+       0-black, 1-r, 2-g,3-b
%                              4-cyan, 5-yellow, 6-purple, 7-white
%

    have_boxlinecolor=~isempty(find(strcmp(textfields,'boxlinecolor')==1));
    if(have_boxlinecolor && isempty(tdata.texts(itxt).boxlinecolor))
        have_boxlinecolor=0;
    end
    if(~have_boxlinecolor)
        text_boxlinecolor=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have boxlinecolor']);
        warning(['Set to default value as ' int2str(text_boxlinecolor)]);
    else
        text_boxlinecolor=tdata.texts(itxt).boxlinecolor;
    end
    fwrite(fid_out,text_boxlinecolor,'int32');


%          +-----------+
%          | INT32     |       Text Box Fill Color
%          +-----------+       0-black, 1-r,2-g,3-b
%                              4-cyan,5-yellow,6-purple,7-white

    have_boxfillcolor=~isempty(find(strcmp(textfields,'boxfillcolor')==1));
    if(have_boxfillcolor && isempty(tdata.texts(itxt).boxfillcolor))
        have_boxfillcolor=0;
    end
    if(~have_boxfillcolor)
        text_boxfillcolor=7;  %default is 7
        warning(['text ' int2str(itxt) ' does not have boxfillcolor']);
        warning(['Set to default value as ' int2str(text_boxfillcolor)]);
    else
        text_boxfillcolor=tdata.texts(itxt).boxfillcolor;
    end
    fwrite(fid_out,text_boxfillcolor,'int32');

%          +-----------+
%          | FLOAT64   |       Angle
%          +-----------+

    have_angle=~isempty(find(strcmp(textfields,'angle')==1));
     if(have_angle && isempty(tdata.texts(itxt).angle))
        have_angle=0;
    end
    if(~have_angle)
        text_angle=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have angle']);
        warning(['Set to default value as ' int2str(text_angle)]);
    else
        text_angle=tdata.texts(itxt).angle;
    end
    fwrite(fid_out,text_angle,'float64');

%          +-----------+
%          | FLOAT64   |       Line Spacing 1-single, 2-double
%          +-----------+                    1.5- one and half
%
    
    have_linespace=~isempty(find(strcmp(textfields,'linespace')==1));
    if(have_linespace && isempty(tdata.texts(itxt).linespace))
        have_linespace=0;
    end
    if(~have_linespace)
        text_linespace=1;  %default is 1
        warning(['text ' int2str(itxt) ' does not have linespace']);
        warning(['Set to default value as ' int2str(text_linespace)]);
    else
        text_linespace=tdata.texts(itxt).linespace;
    end
    fwrite(fid_out,text_linespace,'float64');

%          +-----------+
%          | INT32     |       Text Anchor. 0=left,       1=center,
%          +-----------+                    2=right,      3=midleft   
%                                           4=midcenter,  5=midright, 
%                                           6=headleft    7=headcenter,
%                                           8=headright

    have_anchor=~isempty(find(strcmp(textfields,'anchor')==1));
    if(have_anchor && isempty(tdata.texts(itxt).anchor))
        have_anchor=0;
    end
    if(~have_anchor)
        text_anchor=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have anchor']);
        warning(['Set to default value as ' int2str(text_anchor)]);
    else
        text_anchor=tdata.texts(itxt).anchor;
    end
    fwrite(fid_out,text_anchor,'int32');


%          +-----------+
%          | INT32     |       Zone (0=all)
%          +-----------+

    have_zone=~isempty(find(strcmp(textfields,'zone')==1));
    if(have_zone && isempty(tdata.texts(itxt).zone))
        have_zone=0;
    end
    if(~have_zone)
        text_zone=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have zone']);
        warning(['Set to default value as ' int2str(text_zone)]);
    else
        text_zone=tdata.texts(itxt).zone;
    end
    fwrite(fid_out,text_zone,'int32');


%          +-----------+
%          | INT32     |       Color
%          +-----------+
%                                 0-Black, 1-R, 2-G,3-B
%                                 4-Cyan 5-Yellow
%                                 6-Purple,7-White
%                                 Defalut is 0

    have_color=~isempty(find(strcmp(textfields,'color')==1));
    if(have_color && isempty(tdata.texts(itxt).color))
        have_color=0;
    end
    if(~have_color)
        text_color=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have color']);
        warning(['Set to default value as ' int2str(text_color)]);
    else
        text_color=tdata.texts(itxt).color;
    end
    fwrite(fid_out,text_color,'int32');

%          +-----------+
%          | INT32*N   |       MacroFunctionCommand (string: N = Length + 1)
%          +-----------+
    %this is not supported here. Simply give no commands
    text_MFC_dummy='none';
    plt_write_string(fid_out, text_MFC_dummy);

%          +-----------+
%          | INT32     |       Clipping (0=ClipToAxes,1=ClipToViewport,
%          +-----------+                 2=ClipToFrame)

    have_clipping=~isempty(find(strcmp(textfields,'clipping')==1));
    if(have_clipping && isempty(tdata.texts(itxt).clipping))
        have_clipping=0;
    end
    if(~have_clipping)
        text_clipping=0;  %default is 0
        warning(['text ' int2str(itxt) ' does not have clipping']);
        warning(['Set to default value as ' int2str(text_clipping)]);
    else
        text_clipping=tdata.texts(itxt).clipping;
    end
    fwrite(fid_out,text_clipping,'int32');

%          +-----------+
%          | INT32*N   |       Text.  N=Text Length+1
%          +-----------+
%

    have_str=~isempty(find(strcmp(textfields,'str')==1));
    if(have_str && isempty(tdata.texts(itxt).str))
        have_str=0;
    end
    if(~have_str)
        s=-1;
        warning(['Error : text ' int2str(itxt) ' does not have str']);
        return;
    else
        text_str=tdata.texts(itxt).str;
    end
    plt_write_string(fid_out, text_str);

end  %end if itxt loop

%
%Header section vii. CustomLabel
%

NCustomLabel=0;
have_customlabel=~isempty(find(strcmp(tdatanames,'customlabels')==1)); 
if(have_customlabel && isempty(tdata.customlabbels))
        have_customlabel=0;
end

if(have_customlabel)

   NCustomLabel=0;

   if(iscellstr(tdata.customlabels))
       NCustomLabels=length(tdata.customlabels);
   else
       s=-1;
       display(['Error: customlables must be a cell array (vector) of strings']);
       return;
   end
   if(NCustomLabels>0)
     %Loop through all custom labels if any
%          +-----------+
%          | FLOAT32   |       CustomLabel Marker;  F=599
%          +-----------+
      dummy_float32 = single(599.0);
      fwrite(fid_out,dummy_float32,'float32');

%          +-----------+
%          | INT32     |       Number of labels
%          +-----------+

      fwrite(fid_out,NCustomLabels,'int32');
      for ilabel=1:NCustomLabels
%          +-----------+
%          | INT32*N   |       Text for label 1.  (N=length of label + 1)
%          +-----------+       See note 1.
%          +-----------+
%          | INT32*N   |       Text for label 2.  (N=length of label + 1)
%          +-----------+       See note 1.
%              .
%              .
%              .
%          +-----------+
%          | INT32*N   |       Text for label NumLabels.  
%          +-----------+       (N=length of label + 1) See note 1.

         custom_label=tdata.customlabels{ilabel};
         plt_write_string(fid_out,custom_label);
      end
   end
end

% Header section viii.UserRec

%Loop through all User Records if any

 NUserRecs=0;  %this is hardwired to zero!
               %this code does not accept user defined records
               %this is reserved for future use
 for iusrrec=1:NUserRecs
%          +-----------+
%          | FLOAT32   |       UserRec Marker;  F=699
%          +-----------+
%          +-----------+
%          | INT32*N   |       Text for UserRec.  See note 1.
%          +-----------+
% 
 end

 
% Header section ix. Dataset Auxiliary data.

%Loop through all Dataset auxiliary data (metadata) if any)
NDatasetAux=0;
have_auxdata=~isempty(find(strcmp(tdatanames,'auxdata')==1));
if(have_auxdata && isempty(tdata.auxdata))
   have_auxdata=0;
end

if(have_auxdata)

   NDatasetAux=0;

   if(isstruct(tdata.auxdata))
      NDatasetAux=length(tdata.auxdata);
   else
      s=-1;
      display(['Error: auxdata must be a structure array (vector)']);
      return  
   end

   for iaux=1:NDatasetAux

%          +-----------+
%          | FLOAT32   |       DataSetAux Marker;  F=799.0
%          +-----------+
         dummy_float32 = single(799.0);
         fwrite(fid_out,dummy_float32,'float32');


         %makse sure auxname and auxval are strings
         auxfields=fieldnames(tdata.auxdata(iaux));
         have_auxname=~isempty(find(strcmp(auxfields,'auxname')==1));
         have_auxval= ~isempty(find(strcmp(auxfields,'auxval')==1));

         if(have_auxname && isempty(tdata.auxdata(iaux).auxname))
            have_auxname=0;
         end
         if(have_auxval && isempty(tdata.auxdata(iaux).auxname))
            have_auxval=0;
         end



         if(~have_auxname || ~have_auxval)
            s=-1;
            display(['Error, auxdata ' int2str(iaux) ' must have auxname and auxval']);
            return;
         else
           if(~ischar(tdata.auxdata(iaux).auxname) || ...
              ~ischar(tdata.auxdata(iaux).auxval))
              s=-1;
              display(['Error, auxdata ' int2str(iaux) ' must have auxname and auxval as strings']);
              return;
           else
              aux_name=tdata.auxdata(iaux).auxname;
              aux_val =tdata.auxdata(iaux).auxval;
           end
         end

%          +-----------+
%          | INT32*N   |       Text for Auxiliary "Name".  See note 1.
%          +-----------+

         plt_write_string(fid_out,aux_name);
      
%          +-----------+
%          | INT32     |       Auxiliary Value Format (Currently only
%          +-----------+       allow 0=AuxDataType_String)
          aux_format_dummy=0;
          fwrite(fid_out,aux_format_dummy,'int32');

%          +-----------+
%          | INT32*N   |       Text for Auxiliary "Value".  See note 1.
%          +-----------+

         plt_write_string(fid_out,aux_val);

   end

end


% Header section x. Variable Auxiliary data.

%Loop through all variable's auxdata if any

for ivar=1:tdata.Nvar

    %check if any variable has aux data associated with it
    %if yes, give them, otherwise, skip ??
 

%
%          +-----------+
%          | FLOAT32   |       VarAux Marker;  F=899.0
%          +-----------+
%          +-----------+
%          | INT32*N   |       Variable number (zero based value)
%          +-----------+
%          +-----------+
%          | INT32*N   |       Text for Auxiliary "Name".  See note 1.
%          +-----------+
%          +-----------+
%          | INT32     |       Auxiliary Value Format (Currently only
%          +-----------+       allow 0=AuxDataType_String)
%          +-----------+
%          | INT32*N   |       Text for Auxiliary "Value".  See note 1.
%          +-----------+
%
%
     %this is not yet implemented !


end

%------END of HEADER SECTION-----------------------
%

%This EOH_MARKER separates HEADER SECTION and DATA SECTION in tecplot
%binarly data format
%here EOH means '\eOH' beginning-of-line bindkey
%here EOF means '\eOF' end-of-line bindkey 
%

EOH_MARKER=single(357.0);
fwrite(fid_out, EOH_MARKER ,'float32')

%--------------------------------------------------
% II. DATA SECTION
%
%
%
%  (Repeat the following Sections i, ii, iii for each zone):
%
%######ZONE DATA START########
%   
% Data Section i. for both ordered and Finite Element data
%
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+
%                               Passive Variables
%                               Tecplot can manage many data sets at the same time.
%                               However, within a given data set you must supply the same 
%                               number of variables for each zone. In some cases you may have
%                               data where there are many variables and, for some of the 
%                               zones some of those variables are not important. If that is 
%                               the case, you can set selected variables in those zones to be 
%                               passive. A passive variable is one that will always return 
%                               the value zero if queried (e.g. in a probe) but will not 
%                               involve itself in operations such as the calculations of 
%                               the min and max range. This is very useful when calculating 
%                               default contour levels.
%
%          if "has passive variables" != 0
%               +-----------+
%               | INT32*NV  |  Is variable passive: 0 = no, 1 = yes
%               +-----------+  (Omit entirely if "Has passive variables" 
%                               is 0).
%          +-----------+
%          | INT32     |       Has variable sharing 0 = no, 1 = yes.       
%          +-----------+        
%
%          if "has variable sharing" != 0
%              +-----------+
%              | INT32*NV  |  Zero based zone number to share variable with
%              +-----------+   (relative to this datafile).
%                              (-1 = no sharing).  
%                             (Omit entirely if "Has variable sharing" is 0).
%          +-----------+
%          | INT32     |       Zero based zone number to share 
%          +-----------+       connectivity list with (-1 = no sharing).
%                               FEPOLYGON and FEPOLYHEDRON zones use this 
%                               zone number to share face map data. 
%
%          Compressed list of min/max pairs for each non-shared and 
%          non-passive variable. For each non-shared and non-passive 
%          varaible (as specified above):
%          +-----------+
%          | FLOAT64   |       Min value
%          +-----------+
%          +-----------+
%          | FLOAT64   |       Max value
%          +-----------+
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
%
%
% Data section  ii. specific to ordered zone
%      if "zone number to share connectivity list with" == -1 &&
%         "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined 
%                                   face neighbor connections) * P 
%                              (See note 5 below).
%    
% Data section  iii. specific to Finite element zone
%       if ZoneType is NOT FEPOLYGON or FEPOLYHEDRON:
%          if "zone number to share connectivity lists with" == -1
%                            
%            +-----------+
%            | INT32*N   |   Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%                            This represents JMax sets of adjacency zero
%                            based indices where each set contains L values
%                            where L is:
%                                     2 for LINESEGS
%                                     3 for TRIANGLES
%                                     4 for QUADRILATERALS
%                                     4 for TETRAHEDRONS
%                                     8 for BRICKS
%                            Here JMax means number of Elements
%
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
%       if ZoneType is FEPOLYGON or FEPOLYHEDRON:
%
%           if "zone number to share face map data with" == -1 
%                             
%                             See doc/adkum.pdf page 69 for FaceNumber and
%                             data format regarding faces.
%
%               +-----------+  
%           	| INT32*F   | Face node offsets into the face nodes array 
%               +-----------+ below. Does not exist for FEPOLYGON zones.
%                             F = NumFaces+1. 
%                             This provides starting index in FaceNodes array for
%                             each face.
% 
%             	+-----------+ FaceNodes array
%        		| INT32*FN  | Face nodes array containing the node numbers 
%               +-----------+ for all nodes in all faces. (All faces
%                             must be provided with nodes in order of 
%                             all clockwise or all counter-clockwise order)
%                             FN = total number of face nodes.
%                             
% 
%             	+-----------+ 
%          		| INT32*F   | Elements on the left side of all faces. 
%           	+-----------+ Boundary faces use a negative value which is
%               	      the negated offset into the face boundary 
%                             connection offsets array. A value of "-1"
%               	      indicates there is no left element. 
%                	      F = NumFaces. (give 0 if face is on the edge
%                              of the data)
%                             For 2D data, for each face (edge), place your right-hand 
%                             along the face with your fingers pointing in the
%                             direction of incrementing node numbers 
%                             (i.e. from Node 1 to Node 2). The right side 
%                              of your hand will indicate the right element,
%                              and the left side of your hand will indicate 
%                              the left element
%                              
%                             For 3D data, For each face, curl the fingers of your 
%                             right-hand following the order that the nodes
%                             were presented in the FaceNodes array. Your thumb 
%                             will point to the right element. The left
%                             element is the other adjacent element. If 
%                             the face has more than one neighboring 
%                             element on a single side, you will need to use the 
%                             FaceBoundaryConnectionCounts, FaceBoundaryConnectionElems
%                             and FaceBoundaryConnectionZones array.
%
%                             A negative value in the FaceLeftElems or FaceRightElems 
%                             array is used to indicate that the neighboring
%                             element belongs to another zone. The magnitude of the 
%                             negative number is a pointer to a value in the
%                             FaceBoundaryConnectionElements and 
%                             FaceBoundaryConnectionZones arrays.
%
% 
%      		+-----------+ 
%           | INT32*F   | Elements on the right side of all faces. See 
%      		+-----------+ description of left elements above for more
%               	      details. F = NumFaces.  (give 0 if face is on the edge 
%                             of the data)
%
%
%               if "total number of boundary faces" != 0 
%
%                                   See doc/adkum.pdf page 70 for description
%                                   of Boundary Faces and Boundary Connections
%          		    +-----------+ 
%                   | INT32*NBF |   Boundary face connection offsets into the 
%                   +-----------+   boundary face connecion elements array and
%           				    the boundary face connection zones array.
%            				    The number of elements for a face (F) is
%            				    determined by offset[o] - offset[o-1]
%                       	    where 'o' is the negative value from either
%            				    the left or right elements arrays above. 
%           				    Offset[0] = 0. Offset[1] = 0 so that -1 as
%           				    the left or right element always indicates
%            				    no neighboring element. If the number of 
%           				    elements is 0, then there is no neighboring
%           				    element. 
%            				    NBF = total number of boundary faces + 1. 
% 
%                 +-----------+ 
%                 | INT32*NBI |   Boundary face connection elements. A value of 
%                 +-----------+   "-1" indicates there is no element on part of
%                           	  the face. When working with multiple zones, 
%                                 an additional aspect is folded into the FaceLeftElems and
%                                 FaceRightElems arrays. When the neighboring element is 
%                                 not within the current zone, you cannot identify
%                                 the element by its element number alone. Instead you 
%                                 need to specify both the element number and its
%                                 zone number. This is accomplished using the 
%                                 FaceBoundaryConnectionElements and FaceBoundaryConnectionZones 
%                                 arrays. For each boundary connection, the element 
%                                 number of the boundary connection is stored in the 
%                                 FaceBoundaryConnectionElements array while its zone 
%                                 number is stored in the FaceBoundaryConnectionZones array
%        			  NBI = total number of boundary connections. 
%
%                 +-----------+ 
%          		  | INT32*NBI |   Boundary face connection zones. A value of 
%           	  +-----------+    "-1" indicates the current zone. 
%                                  NBI = total number of boundary connections. 
%
%#########ZONE DATA END###########

%
%find the corresponding matlab format (see "help fread" in matlab)
%for each variable based on vformat
%

have_vformat=~isempty(find(strcmp(tdatanames,'vformat')==1));
if(have_vformat && isempty(tdata.vformat))
   have_vformat=0;
end
%make up empty cell array of formats
var_formatstr=cell([1,tdata.Nvar]);
var_prec=nan*zeros([1,tdata.Nvar]);

for iv=1:tdata.Nvar

   if(have_vformat)
      PrecisoinByCountOfInt32=tdata.vformat(iv)  ;
      if(PrecisoinByCountOfInt32<1 || PrecisoinByCountOfInt32 >6)
         PrecisoinByCountOfInt32=1;  %default to float
         display(['warning: variable number ' int2str(iv) ' does not have proper format']);
         display(['Set to default as 1 (float)']);
      end
   else
      PrecisoinByCountOfInt32=1;
      display(['warning: variable number ' int2str(iv) ' does not have vformat']);
      display(['Set to default as 1 (float)']);
   end

   % PrecisoinByCountOfInt32 is given at beginning of this code
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt  (4 byte)
   %                          4 for shortInt (2 byte)
   %                          5 for Byte
   %                          6 for Bit
   
   var_formatstr{iv}='';   %empty it before assigning values
   switch(PrecisoinByCountOfInt32)
         case 1
           var_formatstr{iv}='float32';  %float,   32 bits (float)
           var_prec(iv)=1;
         case 2
           var_formatstr{iv}='float64';  %float,   64 bits (double)
           var_prec(iv)=2;
         case 3
           var_formatstr{iv}='int32';    %integer, 32 bits (longInt)
           var_prec(iv)=3;
         case 4
           var_formatstr{iv}='int16';    %integer, 16 bits (shortInt)
           var_prec(iv)=4;
         case 5
           var_formatstr{iv}='schar';    %signed char, 8 bits (Byte)
           var_prec(iv)=5;
         case 6
           var_formatstr{iv}='bit1';     %1 bit (Bit)
           var_prec(iv)=6;
   end
end


have_lines=~isempty(find(strcmp(tdatanames,'lines')==1));  %check if have lines
                
%get number of lines in tdata
if(have_lines)  %make sure have lines
    if(isstruct(tdata.lines))  %make sure tdata.lines is structure array
                               %or at least structure 
       Nlines =length(tdata.lines);
    else
       Nlines = 0;
    end
else
    Nlines=0;   
end

%
% Loop through all lines zones
%


%#####################
for iline=1:Nlines

    linefields=fieldnames(tdata.lines(iline)); 
%#########ZONE DATA START###########    
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar
   % var_prec  is calculated above  using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
% 
   %for each variable, output is_passive
   
      % %Complete the following
      
      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %       
      %end
      
end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing: 0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');

if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following
     
     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |       Zero based zone number to share connectivity list with (-1 = no sharing).
%          +-----------+

   %connectivity sharing is not yet implemented
zone_number_to_share_connectivity = -1; % zone number to share connectivity with
                                        % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');



%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+

%--------------% data----------
% write data in the right order
% reshape matrices to 1D vector
%
          %find min and max of each variable

          %check x, y, z, v
          %for lines, first varialble is always x
          %y and z are optional
          %v is used when  tdata.Nvar >=3
          %

           
       have_x=~isempty(find(strcmp(linefields,'x')==1));
       have_y=~isempty(find(strcmp(linefields,'y')==1));
       have_z=~isempty(find(strcmp(linefields,'z')==1));
       have_v=~isempty(find(strcmp(linefields,'v')==1));
       if(have_x && isempty(tdata.lines(iline).x))
          have_x=0;
       end
       if(have_y && isempty(tdata.lines(iline).y))
          have_y=0;
       end
       if(have_z && isempty(tdata.lines(iline).z))
          have_z=0;
       end
       if(have_v && isempty(tdata.lines(iline).v))
          have_v=0;
       end

        if(~have_z)
	      %make sure size of x equal to size of y for (x,y) line (z not given)
	      %set default z value to zero
	      if(~have_x || ~have_y)
		 s=-1;
		 display(['Error: line ' int2str(iline) ' x or y value must be given if you do not give z value']);
		 return;
	      else
		 Ixmax=length(tdata.lines(iline).x);
		 Iymax=length(tdata.lines(iline).y);
		 Imax=min(Ixmax,Iymax); %take minimum value
		 if(Ixmax ~= Iymax)
		     warning(['Warning: line ' int2str(iline) ' x and y length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 x_data=tdata.lines(iline).x(1:Imax);
		 y_data=tdata.lines(iline).y(1:Imax);
		 z_data=zeros(size(x_data)) ; %default to zero for z when it does not exist
	      end
        end
        if(~have_y)
	      %make sure size of x equal to size of z for (x,z) line (y not given)
	      %set default y value to zero
	      if(~have_x || ~have_z)
		 s=-1;
		 display(['Error: line ' int2str(iline) ' x or z value must be given if you do not give y value']);
		 return;
	      else
		 Ixmax=length(tdata.lines(iline).x);
		 Izmax=length(tdata.lines(iline).z);
		 Imax=min(Ixmax,Izmax); %take minimum value
		 if(Ixmax ~= Izmax)
		     warning(['Warning: line ' int2str(iline) ' x and z length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 x_data=tdata.lines(iline).x(1:Imax);
		 z_data=tdata.lines(iline).z(1:Imax);
		 y_data=zeros(size(x_data)) ; %default to zero for y when it does not exist
	      end
        end
        if(~have_x)
	      %make sure size of y equal to size of z for (y,z) line (x not given)
	      %set default x value to zero
	      if(~have_y || ~have_z)
		 s=-1;
		 display(['Error: line ' int2str(iline) ' y or z value must be given if you do not give x value']);
		 return;
	      else
		 Izmax=length(tdata.lines(iline).z);
		 Iymax=length(tdata.lines(iline).y);
		 Imax=min(Izmax,Iymax); %take minimum value
		 if(Izmax ~= Iymax)
		     warning(['Warning: line ' int2str(iline) ' y and z length not equal!']);
		     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
		 end
		 z_data=tdata.lines(iline).z(1:Imax);
		 y_data=tdata.lines(iline).y(1:Imax);
		 x_data=zeros(size(z_data)) ; %default to zero for x when it does not exist
	      end
        end
       
       
        if(have_x && have_y && have_z)
               Ixmax=length(tdata.lines(iline).x);
               Izmax=length(tdata.lines(iline).z);
               Iymax=length(tdata.lines(iline).y);
               Imax=min([Ixmax,Iymax,Izmax]);
               if(Izmax ~= Iymax ||Ixmax ~=Iymax)
                     warning(['Warning: line ' int2str(iline) 'x, y and z length not equal!']);
                     warning(['Only taking the first ' int2str(Imax) ' values, rest of them ignored']);
               end
               x_data=tdata.lines(iline).x(1:Imax);
               y_data=tdata.lines(iline).y(1:Imax);
               z_data=tdata.lines(iline).z(1:Imax);
        end
           %
           %make sure number of variables match == 
           %          tdata.Nvar-number_of_passive_var-number_of_shared_var
           %if not give warning and set to default values
           %
          
           
           if(tdata.Nvar>3)  %look for values in v 
                if(~have_v)
                   s=-1;
                   display(['Error: line ' int2str(iline) ' must have v ']);
                   return
                end
                if(tdata.Nvar>3+size(tdata.lines(iline).v,1)); %number of variables (rows) in v is not enough
                   s=-1;                                       %here give error, 
                                                               %alternatively, we can give warning
                                                               %ans provide zeros for missing variables
                   display(['Error: line ' int2str(iline) ' v must contain at least ' ]);
                   display(['            ' int2str(tdata.Nvar-3) ' variables' ]);
                end
           end 

           %first variable  is always x, if x not given, x is set to zero
           %second variable is always y, if y not given, y is set to zero
           %third variable  is always z, if z not given, z is set to zero
           %rest of from #4 to #tdata.Nvar is obtained from v

           for iv=1:tdata.Nvar
               switch(iv)
                  case 1 %x
                        min_C=min(x_data(:));
                        max_C=max(x_data(:));
                  case 2 %y
                        min_C=min(y_data(:));
                        max_C=max(y_data(:));
                  case 3 %z
                        min_C=min(z_data(:));
                        max_C=max(z_data(:));
                  otherwise %v(iv-3,1:Imax)
                        min_C=min(tdata.lines(iline).v(iv-3,1:Imax));  %fetch the iv-3'th row 'th min
                        max_C=max(tdata.lines(iline).v(iv-3,1:Imax));  %fetch the iv-3'th row 'th max
               end
               fwrite(fid_out,min_C,'float64');
               fwrite(fid_out,max_C,'float64');
           end
%
% zone data
%
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%####################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%####################################################################
%
%

%  write data
  %for ordered nodal, data
  
    for iv=1:tdata.Nvar
        switch(iv)
            case 1 %x
               fwrite(fid_out,x_data,var_formatstr{iv});
            case 2 %y
               fwrite(fid_out,y_data,var_formatstr{iv});
            case 3 %z
               fwrite(fid_out,z_data,var_formatstr{iv});
            otherwise %v(iv-3,1:Imax)
               v_data=tdata.lines(iline).v(iv-3,1:Imax);
               %fwrite(fid_out,tdata.lines(iline).v(iv-3,1:Imax),var_formatstr{iv});
               fwrite(fid_out,v_data,var_formatstr{iv});
        end
    end

%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    % none as misc user def connections = 0 for lines
    
    
% Data section  iii. specific to Finite element zone

%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%#########ZONE DATA END###########
end
%#####################

%
%
% Loop through all surface zones
%
have_surfaces=~isempty(find(strcmp(tdatanames,'surfaces')==1));  %check if have surfaces
%get number of surfaces in tdata
if(have_surfaces)  %make sure have surfaces
    if(isstruct(tdata.surfaces))  %make sure tdata.surfaces is structure array
                                  %or at least structure 
       Nsurfs =length(tdata.surfaces);
    else
       Nsurfs = 0;
    end
else
    Nsurfs=0;   
end

%#####################
for isurf=1:Nsurfs

    surffields=fieldnames(tdata.surfaces(isurf));

%#########ZONE DATA START###########        
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% specify data format for plt file 64 bit floats for each var int32 = 2 per
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar

   % var_prec  is calculated above using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
%
   %for each variable, output is_passive

      % %Complete the following

      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %
      %end

end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing  0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');

if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following

     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |   Zero based zone number to share connectivity 
%          +-----------+   list with (-1 = no sharing).

zone_number_to_share_connectivity = -1; % zone number to share connectivity with
                                        % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');

%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+
      
        %
        %check order (orientation) of surface, x, y, z, v
        %and determin the min and max of each variable 
        %

	%
	%find which order the surface data is organized.
	%
	%If IJ-ordered, then x and y must be
	%given nodal, and have size (Imax*Jmax),Kmax must be one.
	%z will be ignored if not supplied. If z is supplied
	%it will be listed as the third variable

	%If IK-ordered, then x, z must be given as nodal, and have size
	%(IMAX*KMAX), JMAX must be one. y will be ignored if not supplied
	%If y is supplied, it will be listed as the second variable
	%
	%IF JK-ordered, then y,z must be given as nodal, and have size (JMAX*KMAX),
	%IMAX must be one. x will be ignored if not supplied. If x is supplied
	%it will be listed as the first variable

	%check if order is available
	have_order=~isempty(find(strcmp(surffields,'order')==1));
        if(have_order && isempty(tdata.surfaces(isurf).order))
           have_order=0;
        end
	if(have_order)
	    if( isnumeric(tdata.surfaces(isurf).order) ...
	       &&(  tdata.surfaces(isurf).order==1  ... 
		 || tdata.surfaces(isurf).order==2 ...
		 || tdata.surfaces(isurf).order==3) )
		 surface_order=tdata.surfaces(isurf).order;
	    else
		warning(['Oops, surface ' int2str(isurf) ' order value invalid']);
		warning(['set to default value of 3 (IJ-ordered)']);
		surface_order=3;
	    end
	else
	    %default to have order =3 (IJ-ordered data)
        surface_order=3;
        warning(['Oops, surface ' int2str(isurf) ' order not given']);
        warning(['set to default value of 3 (IJ-ordered)']);
	end

	have_x=~isempty(find(strcmp(surffields,'x')==1));
        if(have_x && isempty(tdata.surfaces(isurf).x))
           have_x=0;
        end
	have_y=~isempty(find(strcmp(surffields,'y')==1));
        if(have_y && isempty(tdata.surfaces(isurf).y))
           have_y=0;
        end
	have_z=~isempty(find(strcmp(surffields,'z')==1));
        if(have_z && isempty(tdata.surfaces(isurf).z))
           have_z=0;
        end
	have_v=~isempty(find(strcmp(surffields,'v')==1));
        if(have_v && isempty(tdata.surfaces(isurf).v))
           have_v=0;
        end
	%
	%Check variable size to make sure coordinate variables are nodal and have same size
	%
        switch(surface_order)
	    case 3
		%IJ-ordered, x, y must have same size [Imax, Jmax]
    		if(have_x && have_y)  %must have x and y
        	    Imax_x=size(tdata.surfaces(isurf).x,1);
    		    Jmax_x=size(tdata.surfaces(isurf).x,2);
    		    Imax_y=size(tdata.surfaces(isurf).y,1);
    		    Jmax_y=size(tdata.surfaces(isurf).y,2);
    		    if(Imax_x==Imax_y && Jmax_x==Jmax_y)
                	Imax=Imax_x;
        		Jmax=Jmax_x;
            		Kmax=1;
                else %oops size do not match
                    s=-1;                
                    display(['Error: surface ' int2str(isurf) ' does not have matching x and y dimensions!!']);
                    return;
                end
            else
    		    s=-1;
        	    display(['Error: surface ' int2str(isurf) ' does not have x or y ']);
                return;
            end 
            x_data=tdata.surfaces(isurf).x;
            y_data=tdata.surfaces(isurf).y;
            if(have_z)
               z_data=tdata.surfaces(isurf).z; 
            else
               z_data=zeros(size(x_data));  %z is supplied with zeros
            end 

	    case 2
		%IK-ordered, x, z must have same size [Imax, Kmax]
    		if(have_x && have_z)  %must have x and z
    		    Imax_x=size(tdata.surfaces(isurf).x,1);
    		    Kmax_x=size(tdata.surfaces(isurf).x,2);
    		    Imax_z=size(tdata.surfaces(isurf).z,1);
    		    Kmax_z=size(tdata.surfaces(isurf).z,2);
    		    if(Imax_x==Imax_z && Kmax_x==Kmax_z)
            		Imax=Imax_x;
        			Kmax=Kmax_x;
        			Jmax=1;
    		    else %oops size do not match
        			s=-1;
        			display(['Error: surface ' int2str(isurf) ' does not have matching x and z dimensions!!']);
                	return;
                    end 
                else 
                    s=-1;
                    display(['Error: surface ' int2str(isurf) ' does not have x or z ']);
                    return;
                end 
                x_data=tdata.surfaces(isurf).x;
                z_data=tdata.surfaces(isurf).z;
                if(have_y)
                   y_data=tdata.surfaces(isurf).y;
                else 
                   y_data=zeros(size(x_data));  %y is supplied with zeros
                end 
            
	    case 1
                %JK ordered, y, z must have same size (Jmax, Kmax]
                if(have_y && have_z)  %must have y and z
    		    Jmax_y=size(tdata.surfaces(isurf).y,1);
    		    Kmax_y=size(tdata.surfaces(isurf).y,2);
    		    Jmax_z=size(tdata.surfaces(isurf).z,1);
    		    Kmax_z=size(tdata.surfaces(isurf).z,2);
    		    if(Jmax_y==Jmax_z && Kmax_y==Kmax_z)
        			Jmax=Jmax_y;
        			Kmax=Kmax_y;
        			Imax=1;
    		    else %oops size do not match
        			s=-1;
        			display(['Error: surface ' int2str(isurf) ' does not have matching y and z dimensions!!']);
        			return;
                    end
                else
        	    s=-1;
    		    display(['Error: surface ' int2str(isurf) ' does not have y or z ']);
        	    return;
                end         
                y_data=tdata.surfaces(isurf).y;
                z_data=tdata.surfaces(isurf).z;
                if(have_x)
                   x_data=tdata.surfaces(isurf).x; 
                else
                   x_data=zeros(size(y_data));  % x is supplied with zeros
                end 
        end

           %
           %make sure number of variables match ==
           %          tdata.Nvar-number_of_passive_var-number_of_shared_var
           %if not give warning and set to default values
           %


           if(tdata.Nvar>3)  %look for values in v
                if(~have_v)
                   s=-1;
                   display(['Error: surface ' int2str(isurf) ' must have v ']);
                   return
                end
                if(tdata.Nvar>3+size(tdata.surfaces(isurf).v,1)); %number of variables (rows) in v is not enough
                   s=-1;                                          %here give error,
                                                                  %alternatively, we can give warning
                                                                  %ans provide zeros for missing variables
                   display(['Error: surface ' int2str(isurf) ' v must contain at least ' ]);
                   display(['            ' int2str(tdata.Nvar-3) ' variables' ]);
                end
           end

           %first variable  is always x, if x not given, x is set to zero
           %second variable is always y, if y not given, y is set to zero
           %third variable  is always z, if z not given, z is set to zero
           %rest of from #4 to #tdata.Nvar is obtained from v

           for iv=1:tdata.Nvar
               switch(iv)
                  case 1 %x
                        min_C=min(x_data(:));
                        max_C=max(x_data(:));
                  case 2 %y
                        min_C=min(y_data(:));
                        max_C=max(y_data(:));
                  case 3 %z
                        min_C=min(z_data(:));
                        max_C=max(z_data(:));
                  otherwise %v(iv-3,:,:) %Note that v should be a 3D array
                                         %yet it is addressed here as 2D, in order to get the 
                                         %overall min, max of the 2nd and 3rd dimension
                        min_C=min(tdata.surfaces(isurf).v(iv-3,:));  %fetch the iv-3'th row 'th min
                        max_C=max(tdata.surfaces(isurf).v(iv-3,:));  %fetch the iv-3'th row 'th max
               end
               fwrite(fid_out,min_C,'float64');
               fwrite(fid_out,max_C,'float64');
           end

%
% zone data
%
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%###################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%###################################################################
%
%
%for surfaces, zone_type is 0, i.e. ordered
zone_type = 0;  %ordered

% write data based on datapacking and variable location

	%check if have datapacking
% 	have_datapacking=~isempty(find(strcmp(surffields,'datapacking')==1));
%       if(have_datapacking && isempty(tdata.surfaces(isurf).datapacking))
%          have_datapacking=0;
%       end
%
% 	if(~have_datapacking) %give default value as block (0)
% 	   data_packing=0;
% 	else
% 	   if(   isinteger(tdata.surfaces(isurf).datapacking) ...
% 	      &&   isfinite(tdata.surfaces(isurf).datapacking) ...
% 	      &&  (  tdata.surfaces(isurf).datapacking ==0 ...
% 		      || tdata.surfaces(isurf).datapacking ==1))
% 		    data_packing=tdata.surfaces(isurf).datapacking;
% 	   else
% 	       warning(['datapacking of surface ' int2str(isurf) ' is neither 0 (block) nor 1 (point)!!!']);
% 	       data_packing=0;
% 	       warning(['set default value as 0 (block)']);
% 	   end
%     end
       data_packing=0;

       IsBlock=~data_packing;  %IsBlock=1 (true)  when data_packing=0
                               %IsBlock=0 (false) when data_packing=1
        
	%
	%Wen Long: note that when data_packing=1, i.e. point (IsBlock=0), var
	%location must be nodal!! That is Cell Centered variables can only occur
	%when data_packing is block. Cell Centered variables CAN'T happen when
	%data_packing method is point. For Nodal variables, they can occur for both
	%point packing and block packing methods.
	%

	%check if varloc is available

        have_varloc=~isempty(find(strcmp(surffields,'varloc')==1));
        if(have_varloc && isempty(tdata.surfaces(isurf).varloc))
           have_varloc=0;
        end
        if(~have_varloc) %give default value as 0 (nodal)
           warning(['var location of surface ' int2str(isurf) ' is not given']);
           warning(['set default value as 0 (nodal)']);
           var_loc=0;
           var_specified=0;
        else
           if(     isfinite(tdata.surfaces(isurf).varloc) ...
              &&  (  tdata.surfaces(isurf).varloc ==0 ...
                  || tdata.surfaces(isurf).varloc ==1))
		var_loc=tdata.surfaces(isurf).varloc;
		var_specified=1;

		%Wen Long, make sure when var_specified==1
		%data_packing is 0 (block)
		%rule out conflict conditions (data_packing=1 (point) and var_loc=
		%1 (cell-center) cannot co-exist)
		%That is to say when var_loc=1, data_packing must be zero (block)
		%

		if(data_packing==1 && var_loc==1)
		    s=-1;
		    display(['Error: datapacking of surface ' int2str(isurf) ' (1 point) conflicts with varloc']);
		    return
		end
	   else
	       warning(['var location of surface ' int2str(isurf) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
	       var_loc=0;
	       var_specified=0;
               warning(['set default value as 0 (nodal)']);
           end
        end

    
        %calculate variable location (nodal (0) or cell-centered) for each variable
        surface_varloc=zeros(tdata.Nvar);
        for iv=1:tdata.Nvar
            surface_varloc(iv)=var_loc;
            %make sure for coordinate variables, they are nodal
            switch(surface_order)
               case 3
                    if(iv==1||iv==2)
                       surface_varloc(iv)=0;   %nodal for x and y
                    end
               case 2
                    if(iv==1||iv==3)
                       surface_varloc(iv)=0;   %nodal for x and z
                    end
               case 1
                    if(iv==2||iv==3)
                       surface_varloc(iv)=0;   %nodal for y,and z
                    end
            end
        end

        switch(IsBlock)

            case 0  %point by point 
                    %All variables must be provided as nodal
                    for k=1:Kmax
                    for j=1:Jmax
                    for i=1:Imax
                        for iv=1:tdata.Nvar
                            switch(iv)
                                case 1  %fetch x value
                                   switch(surface_order)
                                        case 3
                                           surface_data_value=x_data(i,j);  %IJ-ordered
                                        case 2
                                           surface_data_value=x_data(j,k);  %JK-orddred
                                        case 1
                                           surface_data_value=x_data(i,k);  %IK-ordered
                                    end
                                case 2  %fetch y value
                                    switch(surface_order)
                                        case 3
                                           surface_data_value=y_data(i,j);  %IJ-ordered
                                        case 2
                                           surface_data_value=y_data(j,k);  %JK-orddred
                                        case 1
                                           surface_data_value=y_data(i,k);  %IK-ordered
                                    end
                                case 3  %fetch z value
                                    switch(surface_order)
                                        case 3
                                           surface_data_value=z_data(i,j);  %IJ-ordered
                                        case 2
                                           surface_data_value=z_data(j,k);  %JK-orddred
                                        case 1
                                           surface_data_value=z_data(i,k);  %IK-ordered
                                    end
                                otherwise %fetch value from v(iv-3,i,j) for IJ-ordered
                                          %                 v(iv-3,j,k) for JK-ordered
                                          %                 v(iv-3,i,k) for IK-ordered
                                    switch(surface_order)
                                        case 3
                                           surface_data_value=tdata.surfaces(isurf).v(iv-3,i,j);  %IJ-ordered
                                        case 2
                                           surface_data_value=tdata.surfaces(isurf).v(iv-3,j,k);  %JK-ordered
                                        case 1
                                           surface_data_value=tdata.surfaces(isurf).v(iv-3,i,k);  %IK-ordered
                                    end
                            end
                            %write the data value based on its precision
                            fwrite(fid_out,surface_data_value,var_formatstr{iv}); 
                        end
                    end
                    end
                    end
            case 1  %block 
                    for iv=1:tdata.Nvar                        
                        switch(surface_order)
                             case 3  %IJ-ordered (Kmax=1)
                                 switch(surface_varloc(iv))
                                      case 0  %nodal

                                        % for k=1:Kmax
                                        % for j=1:Jmax
                                        % for i=1:Imax
                                        %     switch(iv)
                                        %        case 1     %x(1:Imax,1:Jmax,1:Kmax)
                                        %        case 2     %y(1:Imax,1:Jmax,1:Kmax)
                                        %        case 3     %z(1:Imax,1:Jmax,1:Kmax)
                                        %        otherwise  %v(iv-3,1:Imax,1:Jmax)
                                        %     end
                                        % end 
                                        % end
                                        % end

                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Imax,1:Jmax);  
                                             case 2
                                               surface_data_value=y_data(1:Imax,1:Jmax);
                                             case 3
                                               surface_data_value=z_data(1:Imax,1:Jmax);
                                            otherwise
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Imax,1:Jmax); 
                                        end
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv}); 
                                      case 1  %cell centered
                                         %or k=1:Kmax
                                         %or j=1:Jmax-1
                                         %or i=1:Imax-1
                                         %   switch(iv)
                                         %       case 1 %x  %x(1:Imax-1,1:Jmax-1)
                                         %       case 2 %y  %y(1:Imax-1,1:Jmax-1)
                                         %       case 3 %z  %z(1:Imax-1,1:Jmax-1)
                                         %       otherwise  %v(iv-3,1:Imax-1,1:Jmax-1)
                                         %   end
                                         %   fwrite(***)
                                         %end
                                         %end
                                         %end
                                        
                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Imax-1,1:Jmax-1);
                                             case 2
                                               surface_data_value=y_data(1:Imax-1,1:Jmax-1);
                                             case 3
                                               surface_data_value=z_data(1:Imax-1,1:Jmax-1);
                                             otherwise
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Imax-1,1:Jmax-1);

                                        end
                                        %also need to pad zeros to make sure it (only needed for ordered data)
                                        %occupies the same size as nodal
                                        %variables
                                        surface_data_value=reshape(surface_data_value,[1,(Imax-1)*(Jmax-1)]);
                                        surface_data_value=[surface_data_value zeros([1,Imax*Jmax-(Imax-1)*(Jmax-1)])];
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv});
                                 end
                             case 1 %JK-ordered
                                 switch(surface_varloc(iv))
                                      case 0  %nodal

                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Jmax,1:Kmax);
                                             case 2
                                               surface_data_value=y_data(1:Jmax,1:Kmax);
                                             case 3
                                               surface_data_value=z_data(1:Jmax,1:Kmax);
                                             otherwise
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Jmax,1:Kmax);
                                        end
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv});
                                      case 1  %cell centered
                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Jmax-1,1:Kmax-1);
                                             case 2
                                               surface_data_value=y_data(1:Jmax-1,1:Kmax-1);
                                             case 3
                                               surface_data_value=z_data(1:Jmax-1,1:Kmax-1);
                                            otherwise 
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Jmax-1,1:Kmax-1);
                                        end
                                        %also need to pad zeros to make sure it
                                        %occupies the same size as nodal
                                        %variables
                                        surface_data_value=reshape(surface_data_value,[1,(Jmax-1)*(Kmax-1)]);
                                        surface_data_value=[surface_data_value zeros([1,Jmax*Kmax-(Jmax-1)*(Kmax-1)])];                                               
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv});
                                 end
                             case 2 %IK-ordered
                                switch(surface_varloc(iv))
                                      case 0  %nodal
                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Imax,1:Kmax);
                                             case 2
                                               surface_data_value=y_data(1:Imax,1:Kmax);
                                             case 3
                                               surface_data_value=z_data(1:Imax,1:Kmax);
                                             otherwise
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Imax,1:Kmax);
                                        end
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv});
                                      case 1  %cell centered
                                        switch(iv)  %by default matlab fwrite() writes array A[Imax,Jmax,Kmax]
                                                    %with i=1:Imax, j=1:Jmax,k=1:Kmax and with i varying the fastest,
                                                    %j the second and k the last
                                                    %hence  we do not have to have k,j,i nested loops as above
                                             case 1
                                               surface_data_value=x_data(1:Imax-1,1:Kmax-1);
                                             case 2
                                               surface_data_value=y_data(1:Imax-1,1:Kmax-1);
                                             case 3
                                               surface_data_value=z_data(1:Imax-1,1:Kmax-1);
                                             otherwise
                                               surface_data_value=tdata.surfaces(isurf).v(iv-3,1:Imax-1,1:Kmax-1);

                                        end
                                        %also need to pad zeros to make sure it  (only needed for ordered data)
                                        %occupies the same size as nodal
                                        %variables
                                        surface_data_value=reshape(surface_data_value,[1,(Imax-1)*(Kmax-1)]);
                                        surface_data_value=[surface_data_value zeros([1,Imax*Kmax-(Imax-1)*(Kmax-1)])];                                        
                                        fwrite(fid_out,surface_data_value,var_formatstr{iv});
                                end
                        end
                    end
        end

%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    % none as misc user def connections = 0 for lines

    
    
% Data section  iii. specific to Finite Element zone

%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%#########ZONE DATA END###########


%####################
end

%
% Loop through all cube zones
%
have_cubes=~isempty(find(strcmp(tdatanames,'cubes')==1));  %check if have cubes
%get number of cubes in tdata
if(have_cubes)  %make sure have cubes
    if(isstruct(tdata.cubes))     %make sure tdata.cubes is structure array
                                  %or at least structure 
       Ncubes =length(tdata.cubes);
    else
       Ncubes = 0;
    end
else
    Ncubes=0;   
end

%#####################
for icube=1:Ncubes

    cubefields=fieldnames(tdata.cubes(icube));  %obtain field names
                                                %in cubes(icube)

%#########ZONE DATA START###########        
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% specify data format for plt file 64 bit floats for each var int32 = 2 per
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar

   % var_prec  is calculated above  using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
%
   %for each variable, output is_passive

      % %Complete the following

      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %
      %end

end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing: 0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');
if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following

     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |       Zero based zone number to share connectivity list with (-1 = no sharing).
%          +-----------+


zone_number_to_share_connectivity = -1; % zone number to share connectivity with
                                        % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');

%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+

	%
	%for ordered volume data, we need to make sure x,y,z all exist and
	%have the same size (IMAX,JMAX,KMAX)
	%IF any of IMAX, JMAX, KMAX =1, we have to also make sure that
	%the rest of variables have same size and non-centered
	%

        have_x=~isempty(find(strcmp(cubefields,'x')==1));
        have_y=~isempty(find(strcmp(cubefields,'y')==1));
        have_z=~isempty(find(strcmp(cubefields,'z')==1));
        have_v=~isempty(find(strcmp(cubefields,'v')==1));
        if(have_x && isempty(tdata.cubes(icube).x))
           have_x=0;
        end
        if(have_y && isempty(tdata.cubes(icube).y))
           have_y=0;
        end
        if(have_z && isempty(tdata.cubes(icube).z))
           have_z=0;
        end
        if(have_v && isempty(tdata.cubes(icube).v))
           have_v=0;
        end

        %
        %Check variable size to make sure coordinate variables are nodal and have same size
        %
	if(have_x && have_y && have_z)  %must have x, y, z and their sizes must match
	   Imax_x=size(tdata.cubes(icube).x,1);
	   Jmax_x=size(tdata.cubes(icube).x,2);
	   Kmax_x=size(tdata.cubes(icube).x,3);

	   Imax_y=size(tdata.cubes(icube).y,1);
	   Jmax_y=size(tdata.cubes(icube).y,2);
	   Kmax_y=size(tdata.cubes(icube).y,3);

	   Imax_z=size(tdata.cubes(icube).z,1);
	   Jmax_z=size(tdata.cubes(icube).z,2);
	   Kmax_z=size(tdata.cubes(icube).z,3);

	   if(   Imax_x== Imax_y && Imax_x == Imax_z && Imax_y ==Imax_z  ...
	      && Jmax_x== Jmax_y && Jmax_x == Jmax_z && Jmax_y ==Jmax_z  ...
	      && Kmax_x== Kmax_y && Kmax_x == Kmax_z && Kmax_y ==Kmax_z)
		  %good, sizes match

		  Imax=Imax_x;
		  Jmax=Jmax_x;
		  Kmax=Kmax_x;

	   else
	       s=-1;
	       display(['Error: cube ' int2str(icube) ' does not have matching x, y, z coordinate dimensions!! ']);
	       return;
	   end
	else %give error
	   s=-1;
	   display(['Error: cube ' int2str(icube) ' must have 3D coordinates x, y, z !! ']);
	   return;
	end

        %make sure have v when Nvar >3 and also size of v is large enough
        if(tdata.Nvar>3 && ~have_v)
           s=-1;
           display(['Error: cube ' int2str(icube) ' must exist as number of variables > 3']);
           return;
        else
           if(tdata.Nvar > 3+ size(tdata.cubes(icube).v,1))
             s=-1;
             display(['Error: cube ' int2str(icube) ' number of variables (rows) not enough in v']);
             return
           end
        end

        %find min and max of each variable and then output the min and max
        %pair with 'float64' format
        for iv=1:tdata.Nvar
            switch(iv)
               case 1 %x(1:Imax,1:Jmax,1:Kmax)
                   min_C=min(tdata.cubes(icube).x(:)); 
                   max_C=max(tdata.cubes(icube).x(:));
               case 2 %y(1:Imax,1:Jmax,1:Kmax)
                   min_C=min(tdata.cubes(icube).y(:));
                   max_C=max(tdata.cubes(icube).y(:));
               case 3 %z(1:Imax,1:Jmax,1:Kmax)
                   min_C=min(tdata.cubes(icube).z(:));
                   max_C=max(tdata.cubes(icube).z(:));
               otherwise %v(iv-3,:,:,:)  %2nd, 3rd and 4th dimension of v dependent
                                         %on whether variable iv is cell-centered or nodal
                   min_C=min(tdata.cubes(icube).v(iv-3,:));
                   max_C=max(tdata.cubes(icube).v(iv-3,:));
            end
            fwrite(fid_out,min_C,'float64');
            fwrite(fid_out,max_C,'float64');
        end

% zone data
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%###################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%###################################################################

        %write out cube data based on datapacking and variable locations
        %note that x, y, z must be nodal as they are coordinate variables
 
	% check if there is datapacking                                                                                 
% 	 have_datapacking=~isempty(find(strcmp(cubefields,'datapacking')==1));
%        if(have_datapacking && isempty(tdata.cubes(icube).datapacking))
%           have_datapacking=0;
%        end

% 	 if(~have_datapacking) %give default value as block (0)
% 	    data_packing=0;
% 	 else
% 	    if(   isinteger(tdata.cubes(icube).datapacking) ...
% 	       &&  isfinite(tdata.cubes(icube).datapacking) ...
% 	       &&  (  tdata.cubes(icube).datapacking ==0 ...
% 		   || tdata.cubes(icube).datapacking ==1))
% 		 data_packing=tdata.cubes(icube).datapacking;
% 	    else
% 		warning(['datapacking of cube ' int2str(icube) ' is neither 0 (block) nor 1 (point)!!!']);
% 		data_packing=0;
% 		warning(['set default value as 0 (block)']);
% 	    end
% 	 end
 	 data_packing=0;
     
         IsBlock=~data_packing; %IsBlock=1 when data_packing=0 (block)
                                %       =0 when             =1 (point)
	 
	 %check if have varloc
	 have_varloc=~isempty(find(strcmp(cubefields,'varloc')==1));
         if(have_varloc && isempty(tdata.cubes(icube).varloc))
            have_varloc=0;
         end
	 if(~have_varloc) %give default value as 0 (nodal)
	    var_loc=0;
	    var_specified=0;
     else
         tdata.cubes(icube).varloc
	    if(     isfinite(tdata.cubes(icube).varloc) ...
	       &&  (  tdata.cubes(icube).varloc ==0 ...
		     || tdata.cubes(icube).varloc ==1))
		 var_loc=tdata.cubes(icube).varloc;
		 var_specified=1;
		 
		 %Wen Long, make sure when var_specified==1
		 %data_packing is 0 (block) 
		 %rule out conflict conditions (data_packing=1 (point) and var_loc=
		 %1 (cell-center) cannot co-exist)
		 %That is to say when var_loc=1, data_packing must be zero (block)
		 %
		 if(data_packing==1 && var_loc==1)
		     s=-1;
		     display(['Error: datapacking of cube ' int2str(icube) ...
			      ' (1- point) conflicts with varloc (1-centered)']);
		     return
		 end
		 
		 %check if any dimension size is one, if yes, then
		 %make sure we do not have cell-centered data. Volume cells do not
		 %exist if any dimension size degrades to one
		 %
		   
		 if(Imax_x==1 ||Jmax_x==1||Kmax_x==1)
		     if(var_loc==1)
			    s=-1;
			    display(['Error: one dimension of  cube ' int2str(icube) ...
			  	' is singleton and cannot use cell-centered data']);
			    return;
             end 
		 end
        else
		   warning(['var location of cube ' int2str(icube) ' is neither 0 (nodal) nor 1 (cell-center)!!!'])
		   var_loc=0;
		   var_specified=0;
		   warning(['set default value as 0 (nodal)']);
	    end
	 end
	 
	 NV=tdata.Nvar;
         cube_varloc=zeros(tdata.Nvar);
	 for iv=1:tdata.Nvar   %NV is number of variables    
        if(have_varloc)
            if(tdata.cubes(icube).varloc~=1)
    		    cube_varloc(iv) = 0;  %nodal 
        	else    %variable iv is cell centered
                cube_varloc(iv) = 1; %cell centered
            end
        else
            cube_varloc(iv)=0;  %default to nodal
        end
		%make coordinate variables are given nodal 
		%assuming first variable is x, second is y, third is z
		%essentially, x, y, z coordiante variables in this cube
		%must be provided as nodal values
		%rest of variables can be centered if data_packing method is block
		%
		if(iv==1||iv==2||iv==3)
		   cube_varloc(iv)=0;
		end
	 end

         %write out data based on IsBlock and cube_varloc

         switch(IsBlock)
             case 0 %IJK ordered data provided point by point
                    %all variables must be nodal (1:Imax,1:Jmax,1:Kmax)
                      %check size of v (shoud be 4 dimensional, first dimension as variable)
                      %                  2nd, 3rd and 4th as I-,J-,K- (x-,y-,z-) dimension
                      %
                      if(NV>3)
                        Imax_v=size(tdata.cubes(icube).v(iv-3,i,j,k),2);
                        Jmax_v=size(tdata.cubes(icube).v(iv-3,i,j,k),3);
                        Kmax_v=size(tdata.cubes(icube).v(iv-3,i,j,k),4);
                        if(Imax_v~=Imax || Jmax_v~=Jmax || Kmax_v ~= Kmax)
                           s=-1;
                           display(['Error : cube ' int2str(icube) ' v 2nd, 3rd, 4th dimension']);
                           display(['        must be of the same size as x, y, z']);
                           return;
                        end
                      end
                      for k=1:Kmax
                      for j=1:Jmax
                      for i=1:Imax
                          for iv=1:NV
                              switch(iv)
                                 case 1  %x
                                    cube_data_value=tdata.cubes(icube).x(i,j,k);
                                 case 2  %y
                                    cube_data_value=tdata.cubes(icube).y(i,j,k);
                                 case 3  %z
                                    cube_data_value=tdata.cubes(icube).z(i,j,k);
                                 otherwise %iv-3'th variable in v 
                                    cube_data_value=tdata.cubes(icube).v(iv-3,i,j,k);
                              end
                              fwrite(fid_out,cube_data_value,var_formatstr{iv});
                          end
                      end
                      end
                      end
             case 1 %IJK ordered data with block storage
                    for iv=1:NV
                        switch(cube_varloc(iv))
                              case 0 %nodal
                                switch(iv)
                                     case 1  %x
                                        cube_data_value=tdata.cubes(icube).x(1:Imax,1:Jmax,1:Kmax);
                                     case 2  %y
                                        cube_data_value=tdata.cubes(icube).y(1:Imax,1:Jmax,1:Kmax);
                                     case 3  %z
                                        cube_data_value=tdata.cubes(icube).z(1:Imax,1:Jmax,1:Kmax); 
                                     otherwise
                                        Imax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),2);
                                        Jmax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),3);
                                        Kmax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),4);
                                        if(Imax_v~=Imax || Jmax_v~=Jmax || Kmax_v ~= Kmax)
                                           s=-1;
                                           display(['Error : cube ' int2str(icube) ' v 2nd, 3rd, 4th dimension']);
                                           display(['        must be of the same size as x, y, z']);
                                           return;
                                        end
                                        cube_data_value=tdata.cubes(icube).v(iv-3,1:Imax,1:Jmax,1:Kmax);
                                end
                                fwrite(fid_out,cube_data_value,var_formatstr{iv});
                              case 1 %cell-centered
                                
                                switch(iv)
                                     case 1 %x (this will not happen, as x must be nodal as coordinate variable)
                                        cube_data_value=tdata.cubes(icube).x(1:Imax-1,1:Jmax-1,1:Kmax-1);
                                     case 2 %y (this will not happen, as y must be nodal as coordinate variable)
                                        cube_data_value=tdata.cubes(icube).y(1:Imax-1,1:Jmax-1,1:Kmax-1);
                                     case 3 %z (this will not happen, as z must be nodal as coordinate variable)
                                        cube_data_value=tdata.cubes(icube).z(1:Imax-1,1:Jmax-1,1:Kmax-1);
                                     otherwise
                                        Imax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),2);
                                        Jmax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),3);
                                        Kmax_v=size(tdata.cubes(icube).v(iv-3,:,:,:),4);
                                        if(Imax_v~=Imax-1 || Jmax_v~=Jmax-1 || Kmax_v ~= Kmax-1)
                                           s=-1;
                                           display(['Error : cube ' int2str(icube) ' v 2nd, 3rd, 4th dimension']);
                                           display(['        must be of size(x)-1, i.e. [Imax-1,Jmax-1,Kmax-1]']);
                                           return;
                                        end
                                        cube_data_value=tdata.cubes(icube).v(iv-3,1:Imax-1,1:Jmax-1,1:Kmax-1);
                                end
                                %also need to pad zeros to make sure it (only needed for ordered data)
                                %occupies the same size as nodal
                                %variables
                                cube_data_value=reshape(cube_data_value,[1,(Imax-1)*(Jmax-1)*(Kmax-1)]);
                                cube_data_value=[cube_data_value zeros([1,Imax*Jmax*Kmax-(Imax-1)*(Jmax-1)*(Kmax-1)])];
                                fwrite(fid_out,cube_data_value,var_formatstr{iv});
                        end
                    end
         end       

%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    % none as misc user def connections = 0 for lines

    
    
% Data section  iii. specific to Finite element zone

%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%#########ZONE DATA END###########


%####################
end

%
% Loop through all FEline zones
%
%#####################
for iFEline=1:NFElines


    FElinefields=fieldnames(tdata.FElines(iFEline));
%#########ZONE DATA START###########        
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% specify data format for plt file 64 bit floats for each var int32 = 2 per
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar

   % var_prec  is calculated above  using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
%
   %for each variable, output is_passive

      % %Complete the following

      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %
      %end

end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing: 0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');
if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following

     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |       Zero based zone number to share connectivity list with (-1 = no sharing).
%          +-----------+


zone_number_to_share_connectivity = -1; % zone number to share connectivity with
                                        % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');

%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+

        %find Min, Max value based on line order, variable location and x,y,z,v

	%check if have line order

	have_order=~isempty(find(strcmp(FElinefields,'order')==1));
        if(have_order && isempty(tdata.FElines(iFEline).order))
           have_order=0;
        end
	if(~have_order) %default to 4 (3D line)
	    warning(['FEline ' int2str(iFEline) ' does not have line order specified']);
	    FEline_order=4;
	    warning(['set to default as ' int2str(FEline_order) ' (3D line)']);
	else
	    FEline_order=tdata.FElines(iFEline).order;
	end

	if(FEline_order<0 || FEline_order>4)
	   s=-1;
	   display(['Error: FEline ' int2str(iFEline) ' order is incorrect (must be 0,1,2,3 or 4)']);
	   return;
	end

	%
	%Check x,y,z and probe how the lines are defined
	%
	have_x=~isempty(find(strcmp(FElinefields,'x')==1));
	have_y=~isempty(find(strcmp(FElinefields,'y')==1));
	have_z=~isempty(find(strcmp(FElinefields,'z')==1));
        have_v=~isempty(find(strcmp(FElinefields,'v')==1));
        if(have_x && isempty(tdata.FElines(iFEline).x))
           have_x=0;
        end
        if(have_y && isempty(tdata.FElines(iFEline).y))
           have_y=0;
        end
        if(have_z && isempty(tdata.FElines(iFEline).z))
           have_z=0;
        end
        if(have_v && isempty(tdata.FElines(iFEline).v))
           have_v=0;
        end
	%check variable sizes based on order
        switch(FEline_order)
	    case 0  %1D
    		%make sure x exist, and size >=2
    		if(have_x) 
    		    if(~isvector(tdata.FElines(iFEline).x))  
        		s=-1;
        		display(['Error: FEline ' int2str(iFEline) ' x value must be 1D array']);
            		return
	            else
                	Imax_x=length(tdata.FElines(iFEline).x(:)); 
	                if(Imax_x<2)
         	           s=-1;
                	   display(['Error: FEline ' int2str(iFEline) ' x coordiante must have at least 2 points']);
	                   return
        	        else
        		    NumNodes=Imax_x; %get number of nodes
                	end
	            end
	        else 
	            s=-1;
	            display(['Error: FEline ' int2str(iFEline) ' x must exist and be 1D array']);
	            return
	        end 
 	        x_data=tdata.FElines(iFEline).x;
	        if(~have_y)
	          y_data=zeros(size(tdata.FElines(iFEline).x));  %default y to zero if not existing
	        else
                y_data=tdata.FElines(iFEline).y;
                end 
	        if(~have_z) 
	          z_data=zeros(size(tdata.FElines(iFEline).x));  %default z to zero if not existing
	        else
                  z_data=tdata.FElines(iFEline).z;
                end
			
	    case 1  %2D on X, Y coordinates
		    %make sure x, y exist, size match and >=2
		    if(have_x && have_y)
			DIMnumber_x=length(size(tdata.FElines(iFEline).x));
			DIMnumber_y=length(size(tdata.FElines(iFEline).y));
			%if(DIMnumber_x==1&&DIMnumber_y==1)
			if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).y))
			    Imax_x=length(tdata.FElines(iFEline).x(:)); 
			    Imax_y=length(tdata.FElines(iFEline).y(:)); 
			    if(Imax_x~=Imax_y || Imax_x<2 || Imax_y<=2)
				s=-1;
				display(['Error: FEline ' int2str(iFEline) ' x y coordiante must have at least 2 points and size must match']);
				return
			    else
				NumNodes=Imax_x; %get number of nodes
			    end 
			else 
			    s=-1;
			    display(['Error: FEline ' int2str(iFEline) ' x, y must be 1D arrays']);
			    return
			end 
		    else 
			s=-1;
			display(['Error: FEline ' int2str(iFEline) ' x and y must exist and be 1D array']);
			return
		    end 
		    x_data=tdata.FElines(iFEline).x;
		    y_data=tdata.FElines(iFEline).y;
		    if(~have_z)
		       z_data=zeros(size(tdata.FElines(iFEline).x));  %default z to zero if not existing
		    else
                       z_data=tdata.FElines(iFEline).z;
                    end
			
	    case 2  %2D on X, Z coordinates
			%make sure x, z exist, size match and >=2
	  	    if(have_x && have_z)
			   DIMnumber_x=length(size(tdata.FElines(iFEline).x));
		       DIMnumber_z=length(size(tdata.FElines(iFEline).z));
			%if(DIMnumber_x==1&&DIMnumber_z==1)
			if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).z))
			    Imax_x=length(tdata.FElines(iFEline).x(:)); 
			    Imax_z=length(tdata.FElines(iFEline).z(:)); 
			    if(Imax_x~=Imax_z || Imax_x<2 || Imax_z<2)
				s=-1;
				displat(['Error: FEline ' int2str(iFEline) ' x z coordiante must have at least 2 points and size must match']);
				return
			    else 
				NumNodes=Imax_x; %get number of nodes
			    end 
			else
			    s=-1;
			    display(['Error: FEline ' int2str(iFEline) ' x, z must be 1D arrays']);
			    return
			end 
		    else  
			s=-1;
			display(['Error: FEline ' int2str(iFEline) ' x and z must exist and be 1D array']);
                        return
                    end  
                    x_data=tdata.FElines(iFEline).x;
                    z_data=tdata.FElines(iFEline).z;
                    if(~have_y)
                       y_data=zeros(size(tdata.FElines(iFEline).x));  %default y to zero if not existing
                    else
                       y_data=tdata.FElines(iFEline).y;
                    end

	    case 3  %2D on Y, Z coordiantes
                %make sure y, z exist, size match and >=2
                if(have_z && have_y)
                    DIMnumber_z=length(size(tdata.FElines(iFEline).z));
                    DIMnumber_y=length(size(tdata.FElines(iFEline).y));
                    %if(DIMnumber_z==1&&DIMnumber_y==1)
                    if(isvector(tdata.FElines(iFEline).z) && isvector(tdata.FElines(iFEline).y))
                        Imax_z=length(tdata.FElines(iFEline).z(:)); 
                        Imax_y=length(tdata.FElines(iFEline).y(:)); 
                        if(Imax_z~=Imax_y || Imax_z<2 || Imax_y<2)
                            s=-1;
                            display(['Error: FEline ' int2str(iFEline) ' y z coordiante must have at least 2 points and size must match']);
                            return;
                        else
                            NumNodes=Imax_z; %get number of nodes
                        end
                    else 
                        s=-1;
                        display(['Error: FEline ' int2str(iFEline) ' y, z must be 1D arrays']);
                        return ;
                    end   
                else 
                    s=-1;
                    display(['Error: FEline ' int2str(iFEline) ' y and z must exist and be 1D array']);
                    return 
                end 

                if(~have_x)
                   x_data=zeros(size(tdata.FElines(iFEline).y));  %default x to zero if not existing
                else
                    z_data=tdata.FElines(iFEline).z;
                end
                y_data=tdata.FElines(iFEline).y;
                z_data=tdata.FElines(iFEline).z;

	    case 4  %3D on X, Y, Z coordiantes
            %make sure x, y, z exist, size match and >=2
		    if(have_x && have_y && have_z)
			DIMnumber_x=length(size(tdata.FElines(iFEline).x));
			DIMnumber_y=length(size(tdata.FElines(iFEline).y));
			DIMnumber_z=length(size(tdata.FElines(iFEline).z));
			%if(DIMnumber_x==1&&DIMnumber_y==1&&DIMnumber_z==1)
			if(isvector(tdata.FElines(iFEline).x) && isvector(tdata.FElines(iFEline).y) ...
								 && isvector(tdata.FElines(iFEline).z))
			       Imax_x=length(tdata.FElines(iFEline).x(:)); 
			   Imax_y=length(tdata.FElines(iFEline).y(:)); 
			   Imax_z=length(tdata.FElines(iFEline).z(:));  
			    if(  Imax_x~=Imax_y ||Imax_x~=Imax_z ||Imax_y~=Imax_z  ...
				|| Imax_x<2 || Imax_y< 2||Imax_z< 2 )
				s=-1;
				display(['Error: FEline ' int2str(iFEline) ' x y z coordiante must have at least 2 points and size must match']);
				return
			    else
				NumNodes=Imax_x; %get number of nodes
			    end 
			else 
			    s=-1;
			    display(['Error: FEline ' int2str(iFEline) ' x, y, z must be 1D arrays']);
			    return 
			end 
		    else 
			s=-1;
			display(['Error: FEline ' int2str(iFEline) ' x, y, z must exist and be 1D array']);
			return
		    end 

		    x_data=tdata.FElines(iFEline).x;
		    y_data=tdata.FElines(iFEline).y;
		    z_data=tdata.FElines(iFEline).z;
        end 

	%check if varloc is specified
	have_varloc=~isempty(find(strcmp(FElinefields,'varloc')==1));
        if(have_varloc && isempty(tdata.FElines(iFEline).varloc))
           have_varloc=0;
        end
	if(~have_varloc) %give default value as 0 (nodal)
	   var_loc=0;
	   var_specified=0;
	else
	   if(     isfinite(tdata.FElines(iFEline).varloc) ...
	      &&  (  tdata.FElines(iFEline).varloc ==0 ...
		  || tdata.FElines(iFEline).varloc ==1))
		var_loc=tdata.FElines(iFEline).varloc;
		var_specified=1;
		
		%Wen Long, make sure when var_specified==1
		%data_packing is 0 (block) 
		%rule out conflict conditions (data_packing=1 (point) and var_loc=
		%1 (cell-center) cannot co-exist)
		%That is to say when var_loc=1, data_packing must be zero (block)
		%
		if(data_packing==1 && var_loc==1)
		    s=-1;
		    dispplay(['Error: datapacking of FEline ' int2str(iFEline) ' (1 point) conflicts with varloc']);
		    return
		end
	   else
	       warning(['var location of FEline ' int2str(iFEline) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
	       var_loc=0;
	       var_specified=0;
	       warning(['set default value as 0 (nodal)']);
	   end
	end

	   NV=tdata.Nvar;
           FEline_varloc=zeros(tdata.Nvar);
	   for iv=1:NV   %NV is number of variables 
 	       if(~have_varloc)
	            FEline_varloc(iv) = 0;  %nodal by default
	       else    
                  if(tdata.FElines(iFEline).varloc==1)
                      FEline_varloc(iv) = 1; %cell centered
                  else
                      FEline_varloc(iv) = 0  %nodal
                  end
	       end
	       
	       %make coordinate variables nodal
	       %assuming first variable is x, second is y, third is z
	       switch(FEline_order)
		   case 0 %1D line defined on x
		       if(iv==1)
			   FEline_varloc(iv)=0; %x must be nodal
		       end
		   case 1 %2D line defined on (x,y)
		       if(iv==1||iv==2)  
			   FEline_varloc(iv)=0; %x and y must be nodal
		       end
		   case 2 %2D line defined on (x,z)
		       if(iv==1||iv==3) 
			   FEline_varloc(iv)=0; %x and z must be nodal
		       end
		   case 3 %2D line defined on (y,z)
		       if(iv==2||iv==3)
			   FEline_varloc(iv)=0;%y and z must be nodal
		       end
		   case 4 %3D line defined on (x,y,z)
		       if(iv==1||iv==2||iv==3)
			   FEline_varloc(iv)=0; %x,y,and z must be nodal
		       end
	       end
	   end
	 
           %check and make sure v have enough entries when Nvar>3
           if(tdata.Nvar>3 && ~have_v)
               s=-1;
               display(['Error: FEline ' int2str(iFEline) ' must have v to have number of variables >3']);
               return;
           else
               if(tdata.Nvar > 3+ size(tdata.FElines(iFEline).v,1))
                  s=-1;
                  display(['Error: FEline ' int2str(iFEline) ' does not have enough variables (rows)']);
                  return; 
               end
           end

           for iv=1:NV
               switch(iv)
                  case 1 %x
                      min_C=min(x_data(:));
                      max_C=max(x_data(:));
                  case 2 %y
                      min_C=min(y_data(:));
                      max_C=max(y_data(:));
                  case 3 %z
                      min_C=min(z_data(:));
                      max_C=max(z_data(:));

                  otherwise %v(iv-3,:)
                      min_C=min(tdata.FElines(iFEline).v(iv-3,:));
                      max_C=max(tdata.FElines(iFEline).v(iv-3,:));
               end
               fwrite(fid_out,min_C,'float64');
               fwrite(fid_out,max_C,'float64');
           end

%
% zone data
%
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%###################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%###################################################################
%
zone_type=1; %FELINESEG
% write data based on IsBlock (Block vs Point) and variabel location(nodal vs cell-centered)

    %check if have datapacking
%     have_datapacking=~isempty(find(strcmp(FElinefields,'datapacking')==1));
%     if(have_datapacking && isempty(tdata.FElines(iFEline).datapacking))
%        have_datapacking=0;
%     end

%     if(~have_datapacking) %give default value as block (0)
%        data_packing=0;
%     else
%        if(   isinteger(tdata.FElines(iFEline).datapacking) ...
% 	  &&  isfinite(tdata.FElines(iFEline).datapacking) ...
% 	  &&  (  tdata.FElines(iFEline).datapacking ==0 ...
% 	      || tdata.FElines(iFEline).datapacking ==1))
% 	    data_packing=tdata.FElines(iFEline).datapacking;
%        else
% 	   warning(['datapacking of FEline ' int2str(iFEline) ' is neither 0 (block) nor 1 (point)!!!']);
% 	   data_packing=0;
% 	   warning(['set default value as 0 (block)']);
%        end
%     end
    data_packing=0;
    IsBlock=~data_packing ; %IsBlock=0 if data_packing ==1 (POINT)
                            %        1                   0 (POINT)

      %find number of elements from e2n array of the line
      have_e2n=~isempty(find(strcmp(FElinefields,'e2n')==1));
      if(have_e2n && isempty(tdata.FElines(iFEline).e2n))
         have_e2n=0;
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix have to be provided']);
          return;
      else
         NumElements=size(tdata.FElines(iFEline).e2n,1); %
         NumNodes_per_Element=size(tdata.FElines(iFEline).e2n,2); %
         %check if number of nodes per element is two
         if(NumNodes_per_Element~=2) %make sure two nodes per line
                s=-1;
                display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix must have two columns']);
                return;
         end
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEline ' int2str(iFEline) ' e2n (element to node ' ...
                   'connectivity) matrix must have at least 1 row (element)']);
             return;
         end
      end

      switch(IsBlock)
          case 0 % POINT, all variables must be nodal
                  %check v size if Nvar>3
                  if(tdata.Nvar>3)
                     NumNodes_v=size(tdata.FElines(iFEline).v,2);
                     if(NumNodes_v~=NumNodes)
                        s=-1;
                        display(['Error: FEline ' int2str(iFEline) 'v 2nd dimension is not the same as NumNodes']);
                        return;
                     end
                  end

                  for iN=1:NumNodes
                  for iv=1:tdata.Nvar
                      switch(iv)
                          case 1  %x
                             FEline_data_value=x_data(iN);
                          case 2  %y
                             FEline_data_value=y_data(iN);
                          case 3  %z
                             FEline_data_value=z_data(iN);
                          otherwise  %v(iv-3,iN)
                             FEline_data_value=tdata.FElines(iFEline).v(iv-3,iN);
                      end
                      fwrite(fid_out,FEline_data_value,var_formatstr{iv});
                  end
                  end
          case 1 % BLOCK, write variables one by one based on var location
                  for iv=1:tdata.Nvar
                      switch(FEline_varloc(iv))
                        case 0  %nodal
                             switch(iv)
                                case 1 %x
                                   FEline_data_value=x_data(1:NumNodes);
                                case 2 %y
                                   FEline_data_value=y_data(1:NumNodes);
                                case 3 %z
                                   FEline_data_value=z_data(1:NumNodes);
                                otherwise  %v(iv-3,1:NumNodes)
                                   %check v size
                                   NumNodes_v=size(tdata.FElines(iFEline).v,2);
                                   if(NumNodes_v ~= NumNodes)
                                      s=-1;
                                      display(['Error: FEline ' int2str(iFEline) ' v 2nd dimension is not the same as NumNodes']);
                                      return;
                                   end
                                   FEline_data_value=tdata.FElines(iFEline).v(iv-3,1:NumNodes); 
                             end
                             fwrite(fid_out,FEline_data_value,var_formatstr{iv});
                        case 1  %cell center
                             switch(iv)
                                case 1 %x
                                    FEline_data_value=x_data(1:NumElements);
                                case 2 %y
                                    FEline_data_value=y_data(1:NumElements);
                                case 3 %z
                                    FEline_data_value=z_data(1:NumElements);
                                otherwise %v
                                    %check v size
                                    NumElms_v=size(tdata.FElines(iFEline).v,2);
                                   if(NumElms_v ~= NumElements)
                                      s=-1;
                                      display(['Error: FEline ' int2str(iFEline) ' v 2nd dimension is not the same as NumElements']);
                                      return;
                                   end
                                   FEline_data_value=tdata.FElines(iFEline).v(iv-3,1:NumElements);
                             end 
                             %%also need to pad zeros to make sure it  (only needed for ordered data)
                             %%occupies the same size as nodal
                             %%variables
                             %FEline_data_value=[FEline_data_value zeros([1,NumNodes-NumElements])];
                             fwrite(fid_out,FEline_data_value,var_formatstr{iv});
                      end
                  end
    end

%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    % none as misc user def connections = 0 for lines

    
    
% Data section  iii. specific to Finite element zone
if(zone_type>0)
    if(zone_number_to_share_connectivity==-1)
%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).

      %provide zone connectivity based on e2n
     
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+     JMax sets of adjacency zero based indices
%                              eachset contains L values where L is
%                               2 for LINESEGS 
%                               3 for TRIANGLES
%                               4 for QUADRILATERALS
%                               4 for TETRAHEDRONS
%                               8 for BRICKS
      for iE=1:NumElements  %e2n is 1-based, minus-1 to get to 0-based 
                            %node numbers, each row of e2n corresponds to
                            %one element
          fwrite(fid_out,tdata.FElines(iFEline).e2n(iE,:)-1,'int32');
      end

      raw_1to1_face_supplied=0;
      %provide face neighbor relationship of each element's faces
      if(raw_1to1_face_supplied)  %
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
          
      end
      
      NoOfUserDefinedNeighbourConn=0;
      %provide user defined misc face neighbors
      if(NoOfUserDefinedNeighbourConn ~= 0)
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
          
      end 
    end
end
%#########ZONE DATA END###########
end

%####################

%
% Loop through all FEsurface zones
%

%#####################
for iFEsurf=1:NFEsurfs
    
    FEsurffields=fieldnames(tdata.FEsurfaces(iFEsurf));
    
%#########ZONE DATA START###########        
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% specify data format for plt file 64 bit floats for each var int32 = 2 per
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar

   % var_prec  is calculated above  using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
%
   %for each variable, output is_passive

      % %Complete the following

      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %
      %end

end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing: 0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');
if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following

     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |       Zero based zone number to share connectivity list with (-1 = no sharing).
%          +-----------+

zone_number_to_share_connectivity=-1; % zone number to share connectivity with
                                       % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');

%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+


      %find number of elements from e2n array of the surface
      have_e2n=~isempty(find(strcmp(FEsurffields,'e2n')==1));
      if(have_e2n && isempty(tdata.FEsurfaces(iFEsurf).e2n))
         have_e2n=0;
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                   'connectivity) matrix/structure have to be provided']);
          return;
          
      else
                  
         %determine if e2n is a structure or array
         e2n_is_struct=isstruct(tdata.FEsurfaces(iFEsurf).e2n);
         
         if(~e2n_is_struct)  %not a struct
         
            %if array, then find about its number of columns, if 3--triangle,
            %then it is traingular, if 4 then it is quadrilateral, if else
            %give error
            
            DIMnumber=length(size(tdata.FEsurfaces(iFEsurf).e2n));
            if(DIMnumber~=2) %make sure e2n is 2D array
                s=-1;
                display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                         'connectivity) matrix/structure have to be 2D array']); 
                return;
            else  

               %get number of elements
               NumElements=size(tdata.FEsurfaces(iFEsurf).e2n,1);

               %get number of nodes per element
               NumNodes_per_element=size(tdata.FEsurfaces(iFEsurf).e2n,2);

               switch(NumNodes_per_element)
                   case 3
                       zone_type=2; %    2=FETRIANGLE
                   case 4
                       zone_type=3; %    3=FEQUADRILATERAL
                   otherwise
                       s=-1;
                       display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                                '       connectivity) array must have 3 (triangle) or 4 (quadrilateral) columns']); 
                       return;
               end
            end
         
         else %if is a struct then check size of e2n.nodes to make sure
              %it is at least 3.
                            
               if(isvector(tdata.FEsurfaces(iFEsurf).e2n)) %make sure
                                                           %it is a
                                                           %structure
                                                           %vector
                   %get number of elements
                   NumElements=length(tdata.FEsurfaces(iFEsurf).e2n);
                   
                   %check all elements and make sure each element has at
                   %least 3 nodes
                   for icell=1:NumElements
                       %find field names in e2n. Must have nodes field
                       e2nfields=fieldnames(tdata.FEsurfaces(iFEsurf).e2n); 
                       have_nodes=~isempty(find(strcmp(e2nfields,'nodes')==1));
                       if(have_nodes && isempty(tdata.FEsurfaces(iFEsurf).e2n.nodes))
                          have_ndoes=0;
                       end
                       if(have_nodes)
                           
                           %make sure nodes is 1D vector array and also
                           %number of them is greater than 3 (to be
                           %polygon)
                           
                           if(    isvector(tdata.FEsurfaces(iFEsurf).e2n.nodes) ...
                               && length(tdata.FEsurfaces(iFEsurf).e2n.nodes) >=3)
                               %do nothing
                           else
                               s=-1;
                               display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                                    'connectivity) nodes must be 1D array of at least 3 node numbers ']); 
                               return
                           end
                       
                       else %without nodes 
                           s=-1;
                           display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                                    'connectivity) have no valid nodes ']); 
                           return;
                           
                       end
                   end
                   
                   zone_type=6;  %FEPOLYGON if passed all above checks
                   
               else
                   s=-1;
                   display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                            '       connectivity) must be a structure vector']); 
                   return;
               end
         end
          
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEsurface ' int2str(iFEsurf) ' e2n (element to node ' ...
                      '       connectivity) matrix/struct must have at least 1 (element)']);
             return;
         end
      end


      %check if have datapacking
%       have_datapacking=~isempty(find(strcmp(FEsurffields,'datapacking')==1));
%       if(have_datapacking && isempty(tdata.FEsurfaces(iFEsurf).datapacking))
%          have_datapacking=0;
%       end

%       if(~have_datapacking) %give default value as block (0)
%           data_packing=0;
%       else
%          if(   isinteger(tdata.FEsurfaces(iFEsurf).datapacking) ...
%           &&    isfinite(tdata.FEsurfaces(iFEsurf).datapacking) ...
%           &&            (tdata.FEsurfaces(iFEsurf).datapacking ==0 ...
%                       || tdata.FEsurfaces(iFEsurf).datapacking ==1))
%               data_packing=tdata.FEsurfaces(iFEsurf).datapacking;
%          else
%               warning(['datapacking of FEsurface ' int2str(iFEsurf) ' is neither 0 (block) nor 1 (point)!!!']);
%               data_packing=0;
%               warning(['set default value as 0 (block)']);
%          end
%       end
      data_packing=0;
      IsBlock=~data_packing ; %IsBlock=0 if data_packing=1 (Point)
                              %       =1                =0 (Block)
                              
      %find out the order (orientation/placement) of the FE surface

      have_order=~isempty(find(strcmp(FEsurffields,'order')==1));
      if(have_order && isempty(tdata.FEsurfaces(iFEsurf).order))
         have_order=0;
      end
      if(~have_order) %default to 4 (3D surface)
          warning(['FEsurface ' int2str(iFEsurf) ' does not have surface order specified']);
          FEsurf_order=4;
          warning(['set to default as ' int2str(FEsurf_order) ' (3D surface)']);
      else 
          FEsurf_order=tdata.FEsurfaces(iFEsurf).order;
      end

      if(FEsurf_order<1 || FEsurf_order>4)
         s=-1;
         display(['Error: FEsurface ' int2str(iFEsurf) ' order is incorrect (must be 1,2,3 or 4)']);
         return; 
      end

      %
      %check avaibility of x,y,z (coordinate variables) and size consistency
      %based on order of surface. Also calculate number of nodes based on
      %coordinate variables
      %

      have_x=~isempty(find(strcmp(FEsurffields,'x')==1));
      have_y=~isempty(find(strcmp(FEsurffields,'y')==1));
      have_z=~isempty(find(strcmp(FEsurffields,'z')==1));
      have_v=~isempty(find(strcmp(FEsurffields,'v')==1));
      if(have_x && isempty(tdata.FEsurfaces(iFEsurf).x))
         have_x=0;
      end
      if(have_y && isempty(tdata.FEsurfaces(iFEsurf).y))
         have_y=0;
      end
      if(have_z && isempty(tdata.FEsurfaces(iFEsurf).z))
         have_z=0;
      end
      if(have_v && isempty(tdata.FEsurfaces(iFEsurf).v))
         have_v=0;
      end

      switch(FEsurf_order)    
         case 3  %surface is defined on (x,y), e.g. z=f(x,y), v=v(x,y)
            %if z does not exist, it is set to zero
            %(x,y) are coordinate variables, must be nodal
            
            if(have_x && have_y) 
                %make sure x and y are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                          length(tdata.FEsurfaces(iFEsurf).y))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' x and y must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x and y provided on all nodes ']);
               return;
            end
            
            x_data=tdata.FEsurfaces(iFEsurf).x;
            y_data=tdata.FEsurfaces(iFEsurf).y;
            if(~have_z)  %default z value is zero
                z_data=zeros(size(tdata.FEsurfaces(iFEsurf).x));
            else
                z_data=tdata.FEsurfaces(iFEsurf).z;
            end
        case 2 
            %surface is defined on (x,z), e.g. y=f(x,z),v=v(x,z)
            %if y does not exist, it is set to zero
            %(x,z) are coordinate variables, must be nodal
            if(have_x && have_z)                
                %make sure x and z are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                          length(tdata.FEsurfaces(iFEsurf).z))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf)  ...
                        ' x and z must be 1D vectors and have same size']);
                    return;
                end
            else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x and z provided on all nodes ']);
               return;
            end
            x_data=tdata.FEsurfaces(iFEsurf).x;
            if(~have_y)  %default y value is zero
                y_data=zeros(size(tdata.FEsurfaces(iFEsurf).x));
            else
                y_data=tdata.FEsurfaces(iFEsurf).y;
            end
            z_data=tdata.FEsurfaces(iFEsurf).z;

        case 1 %surface is defined on (y,z), e.g. x=f(y,z),v=v(y,z)
            %if x does not exist, it is set to zero
            %(y,z) are coordinate variables, must be nodal
            
             if(have_z && have_y)
                %make sure x and y are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && length(tdata.FEsurfaces(iFEsurf).z)== ...
                           length(tdata.FEsurfaces(iFEsurf).y))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).z);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf) ...
                        ' y and z must be 1D vectors and have same size']);
                    return;
                end
             else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' y and z provided on all nodes ']);
               return;
             end
             z_data=tdata.FEsurfaces(iFEsurf).z;
             y_data=tdata.FEsurfaces(iFEsurf).y;
             if(~have_x)  %default x value is zero
                x_data=zeros(size(tdata.FEsurfaces(iFEsurf).y));
             else
                x_data=tdata.FEsurfaces(iFEsurf).x;
             end
        case 4  %surface is defined on 3D curving surface (x,y,z)
                %v=v(x,y,z)
                %(x,y,z) are coordinate vairables, must be nodal
             if(have_x && have_y && have_z) 
                %make sure x, y and z are vectors and have same size
                if(        isvector(tdata.FEsurfaces(iFEsurf).x)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).y)  ...
                        && isvector(tdata.FEsurfaces(iFEsurf).z)  ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                           length(tdata.FEsurfaces(iFEsurf).y) ...
                        && length(tdata.FEsurfaces(iFEsurf).x)== ...
                           length(tdata.FEsurfaces(iFEsurf).z))
                      NumNodes=length(tdata.FEsurfaces(iFEsurf).x);
                else
                    s=-1;
                    display(['Error: FEsurface ' int2str(iFEsurf)  ...
                        ' x, y and z must be 1D vectors and have same size']);
                    return;
                end
             else
               s=-1;
               display(['Error: FEsurface ' int2str(iFEsurf) ' must have ' ...
                        ' x, y and z provided on all nodes ']);
               return;
             end
             x_data=tdata.FEsurfaces(iFEsurf).x;
             y_data=tdata.FEsurfaces(iFEsurf).y;
             z_data=tdata.FEsurfaces(iFEsurf).z;
      end

      %
      %find variable location using varloc and also order of the surface
      %

      %check if varloc is specified
      have_varloc=~isempty(find(strcmp(FEsurffields,'varloc')==1));
      if(have_varloc && isempty(tdata.FEsurfaces(iFEsurf).varloc))
         have_varloc=0;
      end
      if(~have_varloc) %give default value as 0 (nodal)
        var_loc=0;
        var_specified=0;
      else 
          if(     isfinite(tdata.FEsurfaces(iFEsurf).varloc) ...
                      &&  (tdata.FEsurfaces(iFEsurf).varloc ==0 ...
                        || tdata.FEsurfaces(iFEsurf).varloc ==1))

             var_loc=tdata.FEsurfaces(iFEsurf).varloc;
             var_specified=1;
        
             %Wen Long, make sure when var_specified==1
             %data_packing is 0 (block) 
             %rule out conflict conditions (data_packing=1 (point) and var_loc=
             %1 (cell-center) cannot co-exist)
             %That is to say when var_loc=1, data_packing must be zero (block)
             %
             if(data_packing==1 && var_loc==1)
                s=-1;
                dispplay(['Error: datapacking of FEsurface ' int2str(iFEsurf) ' (1 point) conflicts with varloc']);
              return
             end
          else 
             warning(['var location of FEsurface ' int2str(iFEsurf) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
             var_loc=0;
             var_specified=0;
             warning(['set default value as 0 (nodal)']);
          end 
      end 

      %
      %check and make sure v is available and has enough rows (variables)
      %

      if(tdata.Nvar>3 && ~have_v)
         s=-1;
         display(['Error: FEsurface ' int2str(iFEsurf) ' must have v to have Nvar >3']);
         return;
      else
         if(tdata.Nvar>3)
             if(tdata.Nvar>3+size(tdata.FEsurfaces(iFEsurf).v,1))
                s=-1;
                display(['Error: FEsurface ' int2str(iFEsurf) ' v does not have enough variables(rows)']);
                return;
             end
         end 
      end 
   
      %find variable location for each variable   
      NV=tdata.Nvar;
      FEsurface_varloc=zeros([tdata.Nvar,1]);
      
      for iv=1:NV   %NV is number of variables 

          if(have_varloc)
	      if(tdata.FEsurfaces(iFEsurf).varloc~=1)
	          FEsurface_varloc(iv) = 0;  %nodal 
	      else    %variable iv is cell centered
	          FEsurface_varloc(iv) = 1; %cell centered
	      end 
          else
             FEsurface_varloc(iv)=0;  %default to be nodal
          end
       
          %make coordinate variables are nodal
          %assuming first variable is x, second is y, third is z
          %
          switch(FEsurf_order)
                     
           case 3 %2D surface defined on (x,y)
               if(iv==1||iv==2)  
                   FEsurface_varloc(iv)=0; %x and y must be nodal
               end
               
           case 2 %2D surface defined on (x,z)
               if(iv==1||iv==3) 
                   FEsurface_varloc(iv)=0; %x and z must be nodal
               end
               
           case 1 %2D surface defined on (y,z)
               if(iv==2||iv==3)
                   FEsurface_varloc(iv)=0; %y and z must be nodal
               end
           case 4 %3D surface defined on (x,y,z)
               if(iv==1||iv==2||iv==3)
                   FEsurface_varloc(iv)=0; %x,y,and z must be nodal
               end
          end 
      end 

     %find min and max pair of each variable and output them
     for iv=1:tdata.Nvar
         switch(iv)
           case 1  %x
              min_C=min(x_data(:));
              max_C=max(x_data(:));
           case 2  %y
              min_C=min(y_data(:));
              max_C=max(y_data(:));
            
           case 3  %z
              min_C=min(z_data(:));
              max_C=max(z_data(:));
           otherwise  %v(iv-3,:)
              min_C=min(tdata.FEsurfaces(iFEsurf).v(iv-3,:));
              max_C=max(tdata.FEsurfaces(iFEsurf).v(iv-3,:));
         end 

         fwrite(fid_out,min_C,'float64');
         fwrite(fid_out,max_C,'float64');
     end 

%
% zone data
%
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format
%          +-----------+       as specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%###################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%###################################################################
%
%


switch(IsBlock)
          case 0 % POINT, all variables must be nodal
                  %check v size if Nvar>3
                  if(tdata.Nvar>3)
                     NumNodes_v=size(tdata.FEsurfaces(iFEsurf).v,2);
                     if(NumNodes_v~=NumNodes)
                        s=-1;
                        display(['Error: FEsurface ' int2str(iFEsurf) 'v 2nd dimension is not the same as NumNodes']);
                        return;
                     end
                  end

                  for iN=1:NumNodes
                  for iv=1:tdata.Nvar
                      switch(iv)
                          case 1  %x  %fetch x
                                 FEsurface_data_value=x_data(iN);
                          case 2  %y %fetch y
                                 FEsurface_data_value=y_data(iN);
                          case 3  %z
                                 FEsurface_data_value=z_data(iN);
                          otherwise  %v(iv-3,iN)
                                 FEsurface_data_value=tdata.FEsurfaces(iFEsurf).v(iv-3,iN);
                      end
                      fwrite(fid_out,FEsurface_data_value,var_formatstr{iv});
                  end
                  end
          case 1 % BLOCK, write variables one by one based on var location
                  for iv=1:tdata.Nvar
                        switch(FEsurface_varloc(iv))

                            case 0  %nodal
                                 switch(iv)
                                    case 1 %x
                                       FEsurface_data_value=x_data(1:NumNodes);
                                    case 2 %y
                                       FEsurface_data_value=y_data(1:NumNodes);
                                    case 3 %z
                                       FEsurface_data_value=z_data(1:NumNodes);
                                    otherwise  %v(iv-3,1:NumNodes)
                                       %check v size
                                       NumNodes_v=size(tdata.FEsurfaces(iFEsurf).v,2);
                                       if(NumNodes_v ~= NumNodes)
                                          s=-1;
                                          display(['Error: FEsurface ' int2str(iFEsurf) ' v 2nd dimension is not the same as NumNodes']);
                                          return;
                                       end
                                       FEsurface_data_value=tdata.FEsurfaces(iFEsurf).v(iv-3,1:NumNodes);
                                 end
                                 fwrite(fid_out,FEsurface_data_value,var_formatstr{iv});
                            case 1  %cell center
                                 switch(iv)
                                    case 1 %x
                                       FEsurface_data_value=x_data(1:NumElements);
                                    case 2 %y
                                       FEsurface_data_value=y_data(1:NumElements);
                                    case 3 %z
                                       FEsurface_data_value=z_data(1:NumElements);
                                    otherwise %v
                                        %check v size
                                       NumElms_v=size(tdata.FEsurfaces(iFEsurf).v,2);
                                       if(NumElms_v ~= NumElements)
                                          s=-1;
                                          display(['Error: FEsurface ' int2str(iFEsurf) ' v 2nd dimension is not the same as NumElements']);
                                          return;
                                       end
                                       FEsurface_data_value=tdata.FEsurfaces(iFEsurf).v(iv-3,1:NumElements);
                                 end
                                 %%also need to pad zeros to make sure it  (only needed for ordered data)
                                 %%occupies the same size as nodal
                                 %%variables
                                 %FEsurface_data_value=[FEsurface_data_value zeros([1,NumNodes-NumElements])];
                                 fwrite(fid_out,FEsurface_data_value,var_formatstr{iv});
                        end
                  end
end 


%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    
   % misc user def connections not yet implemented
   % none as misc user def connections = 0 
    
% Data section  iii. specific to Finite element zone

if(zone_type>0)
    if(zone_number_to_share_connectivity==-1)
%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
      %provide zone connectivity based on e2n
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+     JMax sets of adjacency zero based indices
%                              each set contains L values where L is
%                               2 for LINESEGS 
%                               3 for TRIANGLES
%                               4 for QUADRILATERALS
%                               4 for TETRAHEDRONS
%                               8 for BRICKS
%                               >=3 and user defined for FEPOLYGON 
%                               >=4 and user defined for FEPOLYHEDRON

      if(zone_type ~= 6 && zone_type ~=7)   %not FEPOLYGON and not FEPOLYHEDRON
	      for iE=1:NumElements  %e2n is 1-based, minus 1 to get to 0-based 
                                %node numbers, each row of e2n corresponds to
                                %one element
	          fwrite(fid_out,tdata.FEsurfaces(iFEsurf).e2n(iE,:)-1,'int32');
	      end
      else
          for iE=1:NumElements
              fwrite(fid_out,tdata.FEsurfaces(iFEsurf).e2n(iE).nodes(:)-1,'int32');  
              %for POLYGON or POLYHEDRON e2n is a structure
          end
      end

      raw_1to1_face_supplied=0;
      %provide face neighbor relationship of each element's faces
      if(raw_1to1_face_supplied)  %
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).

      end
      NoOfUserDefinedNeighbourConn=0;
      %provide user defined misc face neighbors
      if(NoOfUserDefinedNeighbourConn ~= 0)
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
      end
    end
end


%#########ZONE DATA END###########
end
%####################

%
% Loop through all FEvolume zones
%
%#####################
have_FEvolumes=~isempty(find(strcmp(tdatanames,'FEvolumes')==1));  %check if have FEvolume
%get number of FEvolumes in tdata
if(have_FEvolumes)  %make sure have FEvolume
    if(isstruct(tdata.FEvolumes))  %make sure tdata.surfaces is structure array
                                  %or at least structure 
       NFEvolumes =length(tdata.FEvolumes);
    else
       NFEvolumes = 0;
    end
else
    NFEvolumes =0;   
end

for iFEvol=1:NFEvolumes
       FEvolfields=fieldnames(tdata.FEvolumes(iFEvol));
       
%#########ZONE DATA START###########        
%
% Section i
%
% zone marker = 299.0
%          +-----------+
%          | FLOAT32   |       Zone marker  Value = 299.0
%          +-----------+

dummy_float32 = single(299.0);
fwrite(fid_out,dummy_float32,'float32');

% specify data format for plt file 64 bit floats for each var int32 = 2 per
%          +-----------+
%          | INT32*N   |       variable data format, N=Total number of vars
%          +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
% 

for iv=1:tdata.Nvar

   % var_prec  is calculated above  using tdata.vformat
   % i.e. its value should be 1 for float
   %                          2 for double
   %                          3 for longInt
   %                          4 for shortInt
   %                          5 for Byte
   %                          6 for Bit
   fwrite(fid_out,var_prec(iv),'int32');

end

%passive variables
%          +-----------+
%          | INT32     |       Has passive variables: 0 = no, 1 = yes.
%          +-----------+

    %passive varialbles not implemented yet

has_passive = 0;
fwrite(fid_out,has_passive,'int32');

if(has_passive~=0) %     if "has passive variables" != 0

%            +-----------+
%            | INT32*NV  |     Is variable passive: 0 = no, 1 = yes
%            +-----------+     (Omit entirely if "Has passive variables" is 0).
%
   %for each variable, output is_passive

      % %Complete the following

      %for iv=1:tdata.Nvar
      %      is_passive=0;  %Change this value to 1 if yes for variable iv
      %      fwrite(fid_out,is_passive,'int32');
      %
      %end

end

%variable sharing
%          +-----------+
%          | INT32     |       Has variable sharing: 0 = no, 1 = yes.
%          +-----------+

    %variable sharing is not implemented yet

has_var_share = 0;  % no variable sharing
fwrite(fid_out,has_var_share,'int32');
if(has_var_share)  %if "has variable sharing" != 0
%           +-----------+
%           | INT32*NV  |     Zero based zone number to share variable with (relative to this datafile).
%           +-----------+     (-1 = no sharing).   (Omit entirely if "Has variable sharing" is 0).
     % %Complete the following

     %for iv =1:tdata.Nvar
     %    share_zone_number= 0;  %share variable with zone one in this file
     %                      %1;  %for zone 2
     %                      %2;  %for zone 3, ....
     %                      %-1 for no sharing
     %    fwrite(fid_out,share_zone_number ,'int32');
     %end

end

%
% connectivity sharing
%
%          +-----------+
%          | INT32     |       Zero based zone number to share connectivity list with (-1 = no sharing).
%          +-----------+
zone_number_to_share_connectivity=-1; % zone number to share connectivity with
                                      % -1 for not sharing
fwrite(fid_out,zone_number_to_share_connectivity,'int32');

%
% min and max for each non-shared and non-passive variable
%
%    Compressed list of min/max pairs for each non-shared and non-passive variable. For each
%          non-shared and non-passive varaible (as specified above):
%            +-----------+
%            | FLOAT64   |       Min value
%            +-----------+
%            +-----------+
%            | FLOAT64   |       Max value
%            +-----------+

   %
   %find zone_type based on e2n of this volume. If e2n is array of size NEx4
   %then it is  4=FETETRAHEDRON. If e2n is array of size NEx8, then it is
   %5=FEBRICK. Otherwise, it is 7=FEPOLYHEDRON and faces and element to node 
   %connectivity must be provided in e2n as a structure
   %

      %find number of elements from e2n array of the volume
      have_e2n=~isempty(find(strcmp(FEvolfields,'e2n')==1));
      if(have_e2n && isempty (tdata.FEvolumes(iFEvol).e2n))
         have_e2n=0;
      end
      if(~have_e2n)
          s=-1;
          display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                   'connectivity) matrix/structure have to be provided']);
          return;
          
      else
                  
         %determine if e2n is a structure or array
         e2n_is_struct=isstruct(tdata.FEvolumes(iFEvol).e2n);
         if(~e2n_is_struct)  %not a struct
         
            %if array, then find about its number of columns, if 3--triangle,
            %then it is traingular, if 4 then it is quadrilateral, if else
            %give error
            
            DIMnumber=length(size(tdata.FEvolumes(iFEvol).e2n));
            if(DIMnumber~=2) %make sure e2n is 2D array
                s=-1;
                display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                         'connectivity) matrix/structure have to be 2D array']); 
                return;
            else  

               %get number of elements
               NumElements=size(tdata.FEvolumes(iFEvol).e2n,1);

               %get number of nodes per element
               NumNodes_per_element=size(tdata.FEvolumes(iFEvol).e2n,2);

               switch(NumNodes_per_element)
                   case 4
                       zone_type=4; %    4=FETETRAHEDRON
                   case 8
                       zone_type=5; %    5=FEBRICK
                   otherwise
                       s=-1;
                       display(['Error: FEsurface ' int2str(iFEvol) ' e2n (element to node ' ...
                         'connectivity) array must have 4 (tetrahedron) or 8 (brick) columns']); 
                       return;
               end
               
            end
         
         else %if is a struct then check size of e2n.nodes to make sure
              %it is at least 4.
                            
               if(isvector(tdata.FEvolumes(iFEvol).e2n)) %make sure
                                                         %it is a
                                                         %structure
                                                         %vector
                   %get number of elements
                   NumElements=length(tdata.FEvolumes(iFEvol).e2n);
                   
                   %check all elements and make sure each element has at
                   %least 4 nodes 
                   
                   for icell=1:NumElements
                       %find field names in e2n. Must have nodes field
                       e2nfields=fieldnames(tdata.FEvolumes(iFEvol).e2n); 
                       have_nodes=~isempty(find(strcmp(e2nfields,'nodes')==1));
                       if(have_nodes && isempty(tdata.FEvolumes(iFEvol).e2n.nodes))
                          have_nodes=0;
                       end
                       if(have_nodes)
                           %make sure nodes is 1D vector array and also
                           %number of them is greater than 4 (to be
                           %tetrahedron)
                           
                           if(    isvector(tdata.FEvolumes(iFEvol).e2n.nodes)  ...
                               && length(tdata.FEvolumes(iFEvol).e2n.nodes) >=4)
                               %do nothing
                           else
                               s=-1;
                               display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                                    'connectivity) nodes must be 1D array of at least 4 node numbers ']); 
                               return;
                           end
                       
                       else %without nodes for this element
                           s=-1;
                           display(['Error: Element ' int2str(icell) ' of ' ...
                                    'FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                                    'connectivity) have no valid nodes ']); 
                           return;
                           
                       end
                   end
                   
                   zone_type=7;  %FEPOLYHEDRON if passed all above checks
                   
               else
                   s=-1;
                   display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                     'connectivity) must be a structure vector']); 
                   return;
               end
         end
          
         %check if number of elements is >=1
         if(NumElements <1)
             s=-1;
             display(['Error: FEvolume ' int2str(iFEvol) ' e2n (element to node ' ...
                   'connectivity) matrix/struct must have at least 1 (element)']);
             return;
         end
      end

    %check if have datapacking
%     have_datapacking=~isempty(find(strcmp(FEvolfields,'datapacking')==1));
%     if(have_datapacking && isempty(tdata.FEvolumes(iFEvol).datapacking))
%        have_datapacking=0;
%     end

%     if(~have_datapacking) %give default value as block (0)
%         data_packing=0;
%     else
%         if(      isinteger(tdata.FEvolumes(iFEvol).datapacking) ...
%             &&    isfinite(tdata.FEvolumes(iFEvol).datapacking) ...
%             &&          (  tdata.FEvolumes(iFEvol).datapacking ==0 ...
%                         || tdata.FEvolumes(iFEvol).datapacking ==1))
%                 data_packing=tdata.FEvolumes(iFEvol).datapacking;
%         else
%             warning(['datapacking of FEvolume ' int2str(iFEvol) ' is neither 0 (block) nor 1 (point)!!!']);
%             data_packing=0;
%             warning(['set default value as 0 (block)']);
%         end
%     end
    data_packing=0;
    IsBlock=~data_packing;  %IsBlock=1 (true)  when data_packing=0
                            %IsBlock=0 (false) when data_packing=1

  
    %check if varloc is specified
    have_varloc=~isempty(find(strcmp(FEvolfields,'varloc')==1));
    if(have_varloc && isempty(tdata.FEvolumes(iFEvol).varloc))
        have_varloc=0;
    end
    if(~have_varloc) %give default value as 0 (nodal)
        var_loc=0;
        var_specified=0;
    else
        if(      isfinite(tdata.FEvolumes(iFEvol).varloc) ...
            &&  (  tdata.FEvolumes(iFEvol).varloc ==0 ...
                || tdata.FEvolumes(iFEvol).varloc ==1))

                var_loc=tdata.FEvolumes(iFEvol).varloc;
                var_specified=1;
        
                %Wen Long, make sure when var_specified==1
                %data_packing is 0 (block) 
                %rule out conflict conditions (data_packing=1 (point) and var_loc=
                %1 (cell-center) cannot co-exist)
                %That is to say when var_loc=1, data_packing must be zero (block)
                %
                if(data_packing==1 && var_loc==1)
                    s=-1;
                    dispplay(['Error: datapacking of FEvolume ' int2str(iFEvol) ' (1 point) conflicts with varloc']);
                    return
                end
        else
            warning(['var location of FEvolume ' int2str(iFEvol) ' is neither 0 (nodal) nor 1 (cell-center)!!!']);
            var_loc=0;
            var_specified=0;
            warning(['set default value as 0 (nodal)']);
        end
    end 

    %
    %check avaibility of x,y,z (coordinate variables) and size consistency
    %

    have_x=~isempty(find(strcmp(FEvolfields,'x')==1));
    have_y=~isempty(find(strcmp(FEvolfields,'y')==1));
    have_z=~isempty(find(strcmp(FEvolfields,'z')==1));
    have_v=~isempty(find(strcmp(FEvolfields,'v')==1));
    if(have_x && isempty(tdata.FEvolumes(iFEvol).x))
       have_x=0;
    end
    if(have_y && isempty(tdata.FEvolumes(iFEvol).y))
       have_y=0;
    end
    if(have_z && isempty(tdata.FEvolumes(iFEvol).z))
       have_z=0;
    end
    if(have_v && isempty(tdata.FEvolumes(iFEvol).v))
       have_v=0;
    end


    %volume is defined on 3D coordinates (x,y,z)
    %v=v(x,y,z)
    %(x,y,z) are coordinate vairables, must be nodal
    %and exist
            
    if(have_x && have_y && have_z)                
        %make sure x, y and z are vectors and have same size
        if(        isvector(tdata.FEvolumes(iFEvol).x)  ...
                && isvector(tdata.FEvolumes(iFEvol).y)  ...
                && isvector(tdata.FEvolumes(iFEvol).z)  ...
                &&   length(tdata.FEvolumes(iFEvol).x)== ...
                     length(tdata.FEvolumes(iFEvol).y) ...
                &&   length(tdata.FEvolumes(iFEvol).x)== ...
                     length(tdata.FEvolumes(iFEvol).z))
            NumNodes=length(tdata.FEvolumes(iFEvol).x);
        else
            s=-1;
            display(['Error: FEvolume ' int2str(iFEvol) ...
                     ' x, y and z must be 1D vectors and have same size']);
            return;
        end
        x_data=tdata.FEvolumes(iFEvol).x;
        y_data=tdata.FEvolumes(iFEvol).y;
        z_data=tdata.FEvolumes(iFEvol).z;
    else
        s=-1;
        display(['Error: FEvolume ' int2str(iFEvol) ' must have ' ...
                 ' x, y and z provided on all nodes ']);
        return;
    end

    if(tdata.Nvar>3  &&~have_v) 
        s=-1;
        display(['Error: FEvolume ' int2str(iFEvol) ' must have ' ...
                 '       v for tdata.Nvar >3']);
    end
    
    %
    %give location for each variable in this zone
    %
    FEvol_varloc=zeros(tdata.Nvar);
    for iv=1:tdata.Nvar %NV is number of variables
        
        if(have_varloc)
            if(tdata.FEvolumes(iFEvol).varloc~=1)
                FEvol_varloc(iv) = 0;  %nodal 
            else    %variable iv is cell centered
                FEvol_varloc(iv) = 1; %cell centered
            end
        else
            FEvol_varloc(iv)=0;  %default to nodal
        end
        %make sure coordinate variables are nodal
        %assuming first 3 variables are coordinates for this volume!
        if(iv==1||iv==2||iv==3)  %x, y, z must be nodal
            FEvol_varloc(iv) = 0; %Node
        end
    end 

    %find min and max pair of each variable and output them
     for iv=1:tdata.Nvar
         switch(iv)
           case 1  %x
              min_C=min(x_data(:));
              max_C=max(x_data(:));
           case 2  %y
              min_C=min(y_data(:));
              max_C=max(y_data(:));
           case 3  %z
              min_C=min(z_data(:));
              max_C=max(z_data(:));
           otherwise  %v(iv-3,:)
              min_C=min(tdata.FEvolumes(iFEvol).v(iv-3,:));
              max_C=max(tdata.FEvolumes(iFEvol).v(iv-3,:));
         end 
         fwrite(fid_out,min_C,'float64');
         fwrite(fid_out,max_C,'float64');
     end 
    
%
% zone data
%
%          +-----------+
%          | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
%          +-----------+       specified above.
% 
%  And the data are organized as the following based on zone type,
%  variable location, IsBlock
%   Ordered,
%         Nodal
%              IsBlock=0
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is 
%                     for (K=1;K<=KMax;K++)
%                         for (J=1;J<=JMax;J++)
%                             for (I=1;I<=IMax;I++)
%                                 for(Var=1;Var<=NumVars;Var++)
%                                    Data[Var, I, J, K] = value
%
%             IsBlock=1
%                  Number of values =
%                     IMAX*JMAX*KMAX*NumVars
%                     order is
%                       for (Var=1;Var<=NumVars;Var++)
%                           for (K=1;K<=KMax;K++)
%                              for (J=1;J<=JMax;J++)
%                                  for (I=1;I<=IMax;I++)
%                                     Data[I, J, K, Var] = value;
%
%         Centered
%            IsBlock = 1
%                  Number of values = 
%                    (IMAX-1)*(JMAX-1)*(KMAX-1)*NumVars
%                  Order is: 
%                    for (Var=1;Var<=NumVars;Var++)
%                        for (K=1;K<=(KMax-1);K++)
%                            for (J=1;J<=(JMax-1);J++)
%                                for (I=1;I<=(IMax-1);I++)
%                                    Data[I, J, K, Var] = value;
%            IsBlock = 0
%                  This is Not Allowed!
%
%   Finite element:
%
%        Nodal
%            IsBlock = 0
%                 Number of values = NumNodes*NumVars
%                 Order is:
%                   for (N=1;N<=NumNodes;N++)
%                      for (Var=1;Var<=NumVars;Var++)
%                          Data[Var, N] = value;
%       Nodal
%           IsBlock = 1
%                Number of values = NumNodes*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                       for (N=1;N<=NumNodes;N++)
%                          Data[N, Var] = value;
%       Cell Centered
%           IsBlock = 1
%                Number of values = Jmax(i.e. NumElements)*NumVars
%                Order is:
%                   for (Var=1;Var<=NumVars;Var++)
%                      for (E=1;E<=NumElements;E++)
%                         Data[E, Var] = value;
%           IsBlock = 0
%                This is Not Allowed!
%###################################################################
%====for Centered variable location, IsBlock must be 1 !!!!=========
%###################################################################
%
%

switch(IsBlock)
  case 0 % POINT, all variables must be nodal
	  %check v size if Nvar>3
	  if(tdata.Nvar>3)
	     NumNodes_v=size(tdata.FEvolumes(iFEvol).v,2);
	     if(NumNodes_v~=NumNodes)
		s=-1;
		display(['Error: FEvolume ' int2str(iFEvol) 'v 2nd dimension is not the same as NumNodes']);
		return;
	     end
	  end
	  for iN=1:NumNodes
	  for iv=1:tdata.Nvar
	      switch(iv)
		  case 1  %x  %fetch x
			 FEvolume_data_value=x_data(iN);
		  case 2  %y %fetch y
			 FEvolume_data_value=y_data(iN);
		  case 3  %z
			 FEvolume_data_value=z_data(iN);
		  otherwise  %v(iv-3,iN)
			 FEvolume_data_value=tdata.FEvolumes(iFEvol).v(iv-3,iN);
	      end
	      fwrite(fid_out,FEvolume_data_value,var_formatstr{iv});
	  end
	  end
  case 1 % BLOCK, write variables one by one based on var location
	  for iv=1:tdata.Nvar
		switch(FEvol_varloc(iv))
		    case 0  %nodal
			 switch(iv)
			    case 1 %x
			       FEvolume_data_value=x_data(1:NumNodes);
			    case 2 %y
			       FEvolume_data_value=y_data(1:NumNodes);
			    case 3 %z
			       FEvolume_data_value=z_data(1:NumNodes);
			    otherwise  %v(iv-3,1:NumNodes)
			       %check v size
			       NumNodes_v=size(tdata.FEvolumes(iFEvol).v,2);
			       if(NumNodes_v ~= NumNodes)
				  s=-1;
				  display(['Error: FEvolume ' int2str(iFEvol) ' v 2nd dimension is not the same as NumNodes']);
				  return;
			       end
			       FEvolume_data_value=tdata.FEvolumes(iFEvol).v(iv-3,1:NumNodes);
			 end
			 fwrite(fid_out,FEvolume_data_value,var_formatstr{iv});
		    case 1  %cell center
			 switch(iv)
			    case 1 %x
			       FEvolume_data_value=x_data(1:NumElements);
			    case 2 %y
			       FEvolume_data_value=y_data(1:NumElements);
			    case 3 %z
			       FEvolume_data_value=z_data(1:NumElements);
			    otherwise %v
			       %check v size
			       NumElms_v=size(tdata.FEvolumes(iFEvol).v,2);
			       if(NumElms_v ~= NumElements)
                    s=-1;
                    display(['Error: FEvolume ' int2str(iFEvol) ' v 2nd dimension is not the same as NumElements']);
                    return;
			       end
			       FEvolume_data_value=tdata.FEvolumes(iFEvol).v(iv-3,1:NumElements);
             end
			 fwrite(fid_out,FEvolume_data_value,var_formatstr{iv});
		end
      end
end 


%
% Data section  ii. specific to ordered zone
%
%      if "zone number to share connectivity list with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
%
    % misc user def connections not yet implemented
    % none as misc user def connections = 0 

    
    
% Data section  iii. specific to Finite element zone

if(zone_type>0)
    if(zone_number_to_share_connectivity==-1)
%
%   if "zone number to share connectivity lists with" == -1
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+
%          if "zone number to share connectivity lists with" == -1 &&
%             "raw local 1-to-1 face neighbors are supplied"
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).
%          if "zone number to share connectivity lists with" == -1 &&
%             "number of miscellaneous user defined face neighbor connections" != 0
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
%                              (See note 5 below).
      %provide zone connectivity based on e2n
%            +-----------+
%            | INT32*N   |     Zone Connectivity Data N=L*JMax (see note 3 below).
%            +-----------+     JMax sets of adjacency zero based indices
%                              eachset contains L values where L is
%                               2 for LINESEGS 
%                               3 for TRIANGLES
%                               4 for QUADRILATERALS
%                               4 for TETRAHEDRONS
%                               8 for BRICKS
%                               >=3 and user defined for FEPOLYGON 
%                               >=4 and user defined for FEPOLYHEDRON

      if(zone_type ~= 6 && zone_type ~=7)   %not FEPOLYGON and not FEPOLYHEDRON
          for iE=1:NumElements       %e2n is 1-based, minus-1 to get to 0-based 
                                    %node numbers, each row of e2n corresponds to
                                    %one element
              fwrite(fid_out,tdata.FEvolumes(iFEvol).e2n(iE,:)-1,'int32');
%              tdata.FEvolumes(iFEvol).e2n(iE,:)
          end
      else
          for iE=1:NumElements
              fwrite(fid_out,tdata.FEvolumes(iFEsurf).e2n(iE).nodes(:)-1,'int32');  %for  POLYHEDRON
                                                                                    %e2n is a structure
          end
      end

      raw_1to1_face_supplied=0;
      %provide face neighbor relationship of each element's faces
      if(raw_1to1_face_supplied)  %
%            +-----------+
%            | INT32*N   |     Raw local 1-to-1 face neighbor array.
%            +-----------+     N = (NumElements * NumFacesPerElement)
%                              (See note 4 below).

      end
      NoOfUserDefinedNeighbourConn=0;
      %provide user defined misc face neighbors
      if(NoOfUserDefinedNeighbourConn ~= 0)
%            +-----------+
%            | INT32*N   |     Face neighbor connections.
%            +-----------+     N = (number of miscellaneous user defined face neighbor connections) * P 
      end
    end
end

%#########ZONE DATA END###########
%####################
end

% ------------end of file
fclose(fid_out);

s=0;   %success!
return 


%======================================

function l= plt_write_string(fid_out,string)
% return the size_t of int32
char_hold = int32(string);
l_max = max(size(char_hold));
for ii =1:1:l_max
    fwrite(fid_out,char_hold(ii),'int32');
end
dummy_int32 =0;  %add a null character to the end
fwrite(fid_out,dummy_int32,'int32');
l=l_max+1;



function l= plt_write_numeric(fid_out,numeric1Dvector, numeric_type_str)
% for one dim vector only , dummy nul int32 is not appened
l_max = length(numeric1Dvector);
fwrite(fid_out,numeric1Dvector, numeric_type_str);

% for ii =1:1:l_max
%     fwrite(fid_out,numeric1Dvector(ii), numeric_type_str);
% end
l=l_max+1;


function l= plt_write_dummy(fid_out)
% mean empty
dummy_int32 =0;
fwrite(fid_out,dummy_int32,'int32');
l=1;


% 
% NOTES:
% 
% 1.  All character data is represented by INT32 values.
% 
%      Example:  The letter "A" has an ASCII value of 65.  The WORD
%                written to the data file for the letter "A" is then
%                65.
%                In fortran this could be done by doing the following:
% 
%                Integer*32 I
%                .
%                .
%                I = ICHAR('A');
% 
%                WRITE(10) I
% 
% 
%     All character strings are null terminated (i.e. terminated by a zero value)
% 
% 
% 2.  In FE Data I = Number of points, J = Number of elements, and
%     K = Element type where:
% 
%     0 = Triangles;
%     1 = Quadrilaterals;
%     2 = Tetrahedrons.
%     3 = Bricks.
%     4 = LineSeg
% 
% 
% 3.  This represents JMax sets of adjacency zero based indices where each
%     set contains L values where L is
%     2 for LINESEGS
%     3 for TRIANGLES
%     4 for QUADRILATERALS
%     4 for TETRAHEDRONS
%     8 for BRICKS
% 
% 
% 4.  The raw face neighbor array is dimensioned by number of elements for the
%     zone times the number of faces per element where each member of the array
%     holds the zero based element neighbor of that face. A boundary face is one
%     that has no neighboring element and is represented by a -1. Faces should only
%     be neighbors if they logically share nodes and they should be reciprocal.
%
% 
% 
% 5.  FaceNeighbor Mode   # values  Data
%     ---------------------------------------------------------------------
%     LocalOneToOne       3         cz,fz,cz
%     LocalOneToMany      nz+4      cz,fz,oz,nz,cz1,cz2,...,czn
%     GlobalOneToOne      4         cz,fz,ZZ,CZ
%     GlobalOneToMany     2*nz+4    cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
%     
%     Where:
%         cz = cell in current zone (zero based)
%         fz = face of cell in current zone (zero based)
%         oz = face obsuration flag (only applies to one-to-many):
%                0 = face partially obscured
%                1 = face entirely obscured
%         nz = number of cell or zone/cell associations (only applies to one-to-many)
%         ZZ = remote Zone (zero based)
%         CZ = cell in remote zone (zero based)
%     
%     cz,fz combinations must be unique and multiple entries are
%     not allowed. Additionally, Tecplot assumes that with the
%     one-to-one face neighbor modes a supplied cell face is
%     entirely obscured by it's neighbor.  With one-to-many, the
%     obscuration flag must be supplied.
%     
%     Face neighbors that are not supplied are run through
%     Tecplot's auto face neighbor generator (FE only).
% 
% 5.  Cell centered variable (DATA SECTION)
%     To make reading of cell centered binary data efficient, Tecplot stores
%     IMax*JMax*KMax numbers of cell centered values, where IMax, JMax, and KMax
%     represent the number of points in the I, J, and K directions. Therefore
%     extra zero values (ghost values) are written to the data file for the
%     slowest moving indices. For example if your data's IJK dimensions are 2x3x2
%     a cell centered variable will have 1x2x1 (i.e. (I-1)x(J-1)x(K-1))
%     significant values however 2x3x2 values must be written out because it must
%     include the ghost values. Assume that the two significant cell centered
%     values are 1.5 and 12.5. The ghost values will be output with a zero value.
% 
%     So if the zone was dimensioned 2x3x2 it's cell centered variable would be
%     represented as follows:
%       1.5  0.0    12.5  0.0     0.0  0.0     0.0  0.0     0.0  0.0     0.0  0.0
% 
%     If the zone was dimensioned 3x2x2 it's cell centered variable would be
%     represented as follows:
%       1.5  12.5  0.0     0.0  0.0  0.0     0.0  0.0  0.0     0.0  0.0
%       0.0
% 
%     and if the zone was dimensioned 2x2x3 it's cell centered variable would be
%     represented as follows:
%       1.5  0.0    0.0  0.0     12.5  0.0     0.0  0.0     0.0  0.0     0.0  0.0
% 
%     Obviously for large variables the wasted space is less significant that it
%     is for the small example above.
%-----------------------------------------------------------------------

