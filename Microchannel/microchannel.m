%% Image Processing for Microchannel Flow experiment
% Prototype code:
% Using images taken in the experiment:
%
% 1st image fluid flow: 20 mu L/min
% 2st image fluid flow: 10 mu L/min
%
% Microchannel 1: 1808.bmp & 1813.bmp
% Microchannel 2: 1824.bmp & 1828.bmp
% Microchannel 3: 1834.bmp & 1838.bmp
% Microchannel 4: 1843.bmp & 1847.bmp
% Microchannel 5: 1855.bmp & 1860.bmp

target = {'1808.bmp' '1813.bmp'...
    '1824.bmp' '1828.bmp' ...
    '1834.bmp' '1838.bmp' ...
    '1843.bmp' '1847.bmp' ...
    '1855.bmp' '1860.bmp'};
%or 
number = [1808 1813 1824 1828 1834 1838 1843 1847 1855 1860];

%Pix_value = [0.31 0.31 0.31 0.31 0.31 0.31 0.31 0.31 0.31 0.31];

%% Load images
name = strcat(num2str(number(1)),'.bmp');
image = imread(name);

%% Apply filter
h = fspecial('unsharp');
image2 = imfilter(image,h); % apply X filter 

%% Convert to gray scale
image3 = rgb2gray(image); % conver Original to gray scales
image4 = rgb2gray(image2); % conver Filtered to gray scales

%% Create image negative
% Original Gray
image5 = double(image3);
scale_image = 1/max(image5(:)); % scale_image = 1/255
neg_image5 = 1-image5*scale_image;

% Filtered Gray
image6 = double(image4);
scale_image = 1/max(image6(:)); % scale_image = 1/255
neg_image6 = 1-image6*scale_image;

%% Remove values below an specific Pix_value
Pix_value = 0.31;

[n,m,p] = size(neg_image5);
for k = 1:p
    for j = 1:n;
        for i = 1:m;
            if neg_image5(j,i,k) >= Pix_value
                neg_image5(j,i,k) = 0;
            else
                % do nothing
            end
        end
    end
end
[n,m,p] = size(neg_image6);
for k = 1:p
    for j = 1:n;
        for i = 1:m;
            if neg_image6(j,i,k) >= Pix_value
                neg_image6(j,i,k) = 0;
            else
                % do nothing
            end
        end
    end
end

%% Plot pictures
figure(1); 
subplot(2,3,1); imshow(image); title('Original');
subplot(2,3,4); imshow(image2); title('Filtered Image');
subplot(2,3,2); imshow(image3); title('Gray Image with filter');
subplot(2,3,5); imshow(image4); title('Gray Image with out filter');
subplot(2,3,3); imshow(neg_image5); title('Original, filtered an processed');
subplot(2,3,6); imshow(neg_image6); title('negative of the filtered and processed');

%% Save gray_imageX
name2 = strcat(num2str(number(1)),'_gray.bmp');
imwrite(image4,name2,'bmp');

%% Save neg_imageX
name2 = strcat(num2str(number(1)),'_neg.bmp');
imwrite(neg_image6,name2,'bmp');

%% Compute Mixing Index
I = neg_image6;      % Use values in neg_image6 for the difusion analysis.
I_mean = mean(I(:)); % Find the average of the image intensity values.
x = find(I);         % find indices of non-zero values.
N = length(x);       % Total number of non-zero values.

Mix_index = 1-(1/I_mean)*sqrt(sum(sum((I - I_mean).^2))/N);

% Now that I sure what I want to do I'll create a function to evaluate
% automatically data picture.