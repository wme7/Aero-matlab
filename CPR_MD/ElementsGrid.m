%% Build Grid
clc; clear all; close all;

numberSPs = 4;

[xi,w] = GaussLegendre(numberSPs);

normSPs = xi;

         range = [0,1];
numberElements = 9;
   numberNodes = numberElements*numberSPs;
  nodedbCoords = linspace(range(1),range(2),numberElements+1); 
   elementSize = (range(2)-range(1))/numberElements;
 elementCenter = (nodedbCoords(1:numberElements) + ...
                  nodedbCoords(2:numberElements+1)) / 2;
        metric = elementSize/2;
[scaledSPs,xc] = meshgrid(elementCenter,metric*normSPs);
    nodeCoords = scaledSPs+xc; % x-Grid
  elementNodes = reshape(1:numberNodes,numberSPs,numberElements);
elementdbNodes = [elementNodes(1,:); elementNodes(numberSPs,:) ];