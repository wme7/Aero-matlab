function[DepMat RgbMat] = Kinect()
%---------------------------------------------------------------
% This function retrive Dep and Rgb color matries from mex code
% And return as DepMat RgbMat, respectly, after reshape them.
%---------------------------------------------------------------
[a b] = getimagedata(2);
DepMat = reshape(a,640,480)';
temp = (reshape(b(1:3:end),640,480))';
RgbMat(:,:,1) = temp(:,end:-1:1);
temp = (reshape(b(2:3:end),640,480))';
RgbMat(:,:,2) = temp(:,end:-1:1);
temp = (reshape(b(3:3:end),640,480))';
RgbMat(:,:,3) = temp(:,end:-1:1);
end
