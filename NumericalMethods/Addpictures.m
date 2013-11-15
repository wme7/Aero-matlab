%% Load Images
wallpaper1=imread('Nike.jpg');
wallpaper2=imread('Fire.jpg');

%% Display Wallpaper1
image(wallpaper1);
axis equal; grid off;

%% Display Wallpaper2
image(wallpaper2);
axis equal; grid off;

%% Scale and Display Wallpaper1
darkerlogo=0.4*wallpaper1;
image(darkerlogo);
axis equal; grid off;

%% Add images
% Image fusion is made here ;)
fusion=darkerlogo+wallpaper2;
image(fusion);
axis equal; grid off;