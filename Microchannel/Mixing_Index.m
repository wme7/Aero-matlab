% Mixing Index of two fluids in a Microchannel.
% Experiments on Applied Mechanics.
% By Manuel Diaz 2012.10.20.

%% INPUT IMAGES & CORRECTION VALUES
Im_number = [1808 1813 1824 1828 1834 1838 1855 1860 1843 1847]; 
Pix_value = [0.28 0.28 0.295 0.295 0.285 0.285 0.30 0.30 0.30 0.30]; 
% The Pix_value was choosen by observing the result of the filter process
% until we find an optimun condition to visualize the intesity of the
% mixing region. 

N_Im = length(Im_number);
N_Pv = length(Pix_value);
if N_Im ~= N_Pv 
  Error('N_Im and N_Pv have not the same length');
  break
end

%% MAIN LOOP
for ii = 1:N_Pv
    %% Load Image
    name = strcat(num2str(Im_number(ii)),'.bmp');
    image = imread(name);
    
    %% Apply filter and convert to Gray Intensity scale.
    h = fspecial('unsharp');    % load 'x' filter
    image2 = imfilter(image,h); % apply 'x' filter
    image4 = rgb2gray(image2); % conver Filtered to gray scales
    
    %% Create Negative version of the Gray-Filtered image.
    image6 = double(image4);
    scale_image = 1/max(image6(:)); % scale_image = 1/255
    neg_image6 = 1-image6*scale_image;
    
    [n,m,p] = size(neg_image6);
    for k = 1:p
        for j = 1:n;
            for i = 1:m;
                if neg_image6(j,i,k) >= Pix_value(ii)
                    neg_image6(j,i,k) = 0;
                else
                    % do nothing
                end
            end
        end
    end
    
    %% Save gray_imageX
    name2 = strcat(num2str(Im_number(ii)),'_gray.bmp');
    imwrite(image4,name2,'bmp');
    
    %% Save neg_imageX
    name2 = strcat(num2str(Im_number(ii)),'_neg.bmp');
    imwrite(neg_image6,name2,'bmp');
    
    %% Compute Mixing Index
    I = neg_image6;      % Use values in neg_image6 for the difusion analysis.
    I_mean = mean(I(:)); % Find the average of the image intensity values.
    x = find(I);         % find indices of non-zero values.
    N = length(x);       % Total number of non-zero values.
    
    Mix_index(ii) = 1-(1/I_mean)*sqrt(sum(sum((I - I_mean).^2))/N);
    
end
%% Plot results
figure(1)
% Velocity: 20 mu L/min
subplot(1,2,1)
v1 = Mix_index(1:2:9);
bar(abs(v1)); title('Abs of Diffusion coef @ v_1');
set(gca,'XTickLabel',{'a','b','c','d','e'})
ylabel('Mixing Index');

% Velocity: 10 mu L/min
subplot(1,2,2);
v2 = Mix_index(2:2:10);
bar(abs(v2)); title('Abs of Diffusion coef @ v_2');
set(gca,'XTickLabel',{'a','b','c','d','e'})
ylabel('Mixing Index');
