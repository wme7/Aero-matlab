%% How to create an animated gif picture in Matlab:
clear; close all;

Z = peaks; % using built in function peaks
surf(Z)    % creates our 3D surface
axis tight % we ask for a tight space to get a verter view of our result
set(gca,'nextplot','replacechildren','visible','on') % gca -> grid off!

f = getframe; % the take a snapshot of the image
% f is a structure file: f.cdata and f.colormap <- keep this in mind!

% transfom a RGB image into a indexed image using 256 colors
[im,map] = rgb2ind(f.cdata,256,'nodither'); 
% We then get two results:
% 1. "im" is our indexed image/frame with format uint*. 
%         *Notice that is a 3d vector: im(1,1,1)
% 2. "map" is our new map in 256 colors that matches approximately the 
%          original RBG colors. 
% IMPORTANT: use the option 'nodither' otherwize is going to take forever!

n = 20; % -> I define the amount of frames that my image will have. 
% Normally we don't know this in a real numerical experiment.
% But, if we know it, we're ought to set it form the very begining!

% Now this is a nice trick! 
% to avoid Matlab message of growing variable over time we just set 
% the last layer of our frame temporaly to zero.
im(1,1,1,n) = 0;

% Now the loog with the computation of our next result and it's frame:
for k = 1:n
  surf(cos(2*pi*k/20)*Z,Z)
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither'); % <- look closely now!  
  % we just converted the RGB image to an indexed image again... but! 
  % "f.cdata" we use instead our predefined 256 colormap: "map". 
end

% Now we are ready to write our animated gif file:
imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf)
% see matlab reference for more imwrite options.

% Modifications and comments added by Manuel Diaz