 % Copyright 2001, Brown University, Providence, Rhode Island.
 %
 % All Rights Reserved
 % 
 % Permission to use this software for noncommercial research and
 % educational purposes is hereby granted without fee.
 % Redistribution, sale, or incorporation of this software into a
 % commercial product is prohibited.
 % 
 % BROWN UNIVERSITY DISCLAIMS ANY AND ALL WARRANTIES WITH REGARD TO
 % THIS SOFTWARE,INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 % AND FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN
 % UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 % DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE,
 % DATA OR PROFITS.

 % This routine reads in a .neu (neutral Gambit file)
 % it basically requires:
 %                         umNel        = number of elements
 %                         umNnodes     = number of unique vertices
 %                         umVertX      = X coordinates of vertices
 %                         umVertY      = Y coordinates of vertices
 %                         umBoundaryNodeType = integer flag for different
 %                                              types of boundary conditions
 %                                              at node
 %                         elmttonode   = [umNelx3] = list of vertices in each
 %                                                            element
 % It creates:             umElmtToElmt = [umNelx3] = element to element list
 %                         umElmtToFace = [umNelx3] = element to face list
 %
 % It also shuffles the element order so that the graph of the element 
 % connectivities is bandwidth minimized (using Reverse Cuthill McKee)



 umNfaces   = 3;
 umSTOPNOW  = 0;
 umVnum    = [[1,2];[2,3];[1,3]];

% Set file name for input
 if(~exist('umFileName'))
  umFileName = 'block2.neu';
 end
 umFid = fopen(umFileName, 'rt');

 % read intro 
 for i=1:6 
  line = fgetl(umFid);
 end

 % fine number of nodes and number of elements
 dims = fscanf(umFid, '%d');

 umNnodes = dims(1);
 umNel    = dims(2);
 
 for i=1:2 
  line = fgetl(umFid);
 end

 % read node coordinates
 umVertX = (1:umNnodes);
 umVertY = (1:umNnodes);

 for i = 1:umNnodes
  line = fgetl(umFid);
  tmpx = sscanf(line, '%lf');
  umVertX(i) = tmpx(2);
  umVertY(i) = tmpx(3);
 end
  
 for i=1:2 
  line = fgetl(umFid);
 end

 % read element to node connectivity
 elmttonode = zeros(umNel, 3);

 for k = 1:umNel
  line   = fgetl(umFid);
  tmpcon = sscanf(line, '%lf');
  elmttonode(k,1) = tmpcon(4);
  elmttonode(k,2) = tmpcon(5);
  elmttonode(k,3) = tmpcon(6);
 end

