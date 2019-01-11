%% Disparity map computation with local methods (SSD cost)

% The 'stereo_computation' function computes the disparity between a pair
% of rectified images using a local method based on a matching cost 
% between two local windows.
%
% In this part you are asked to implement the SSD cost
%

% Data images (rectified images)
Irgb{1} = imread('Data/scene1.row3.col3.ppm');
Irgb{2} = imread('Data/scene1.row3.col2.ppm');

% Disparity ground truth image
Igt = imread('Data/truedisp.row3.col3.pgm');

% Set minimum and maximum disparity
min_dis=0;
max_dis=16;
% Set window (patch) size
wsize=21;
% Choose cost function
cost_function='SSD';
% Question 1: 
% Complete the 'stereo_computation' function with the
% computation of the SSD cost
disparity = stereo_computation(Irgb{1}, Irgb{2}, min_dis, max_dis, wsize, cost_function,0);
figure;
title('disparity')
imshow(disparity,[min_dis max_dis]);

% Display the error image with respect to the ground truth
[h,w] = size(Igt);
figure;
imagesc(abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis));
colorbar;

% Question 2:
% How do the error and the occluded areas are related?

% Question 3:
% Compute the global error
error = abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis);
[h,w] = size(error);
globalError = sum(sum(error))/(h*w);

% Question 4:
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 21x21,
% 31x31). Comment the results and differences in the results.

%% NCC similarity function

% Question 5:
% Complete the 'stereo_computation' function with the
% computation of the NCC cost.

wsize=3;
min_dis=0;
max_dis=16;
disparity = stereo_computation(Irgb{1}, Irgb{2}, min_dis, max_dis, wsize, 'NCC',0);
figure;
imshow(disparity,[min_dis max_dis]);

[h,w] = size(Igt);
figure;
imagesc(abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis));
colorbar;

error = abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis);
[h,w] = size(error);
globalError = sum(sum(error))/(h*w);                                                %#ok<NASGU>

% Question 6:
% Compare the results with those obtained with the SDD cost.


%% Bilateral weights

% Question 7:
% Complete the 'stereo_computation' function with the
% bilateral weights

% Data images (rectified images)
Irgb{1} = imread('Data/scene1.row3.col3.ppm');
Irgb{2} = imread('Data/scene1.row3.col2.ppm');

% Disparity ground truth image
Igt = imread('Data/truedisp.row3.col3.pgm');

wsize=25;
min_dis=0;
max_dis=16;
[disparity,cost,weight] = stereo_computation(Irgb{1}, Irgb{2}, min_dis, max_dis, wsize, 'SSD',1);
weight ;
figure;
imshow(disparity,[min_dis max_dis]);

[h,w] = size(Igt);
figure;
imagesc(abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis));
colorbar;

error = abs(double(Igt(18:h-18,18:w-18))-disparity(18:h-18,18:w-18)*255/max_dis);
[h,w] = size(error);
globalError = sum(sum(error))/(h*w);

% Question 8:
% Evaluate the results changing the window size (e.g. 5x5, 9x9, 31x31)
% and compare to the previous case that uses uniform weights (SDD cost).

%% Disparity with lab 5 images


% Data images (rectified images)
Irgb{1} = imread('Data/output_0.png');
Irgb{2} = imread('Data/output_1.png');

scale = 0.25;
I1 = imresize(Irgb{1}, scale);
I2 = imresize(Irgb{2}, scale);

min_dis = -100;
max_dis = 100;
[disparity, cost] = stereo_computation(I1, I2, min_dis, max_dis, 15, 'SSD',1);
figure;
imshow(abs(disparity),[0 max_dis]);colorbar;

% Display only foreground objects
Mask = imread('Data/mask.png');
Mask = imresize(Mask, scale); 
Mask = double(255 - Mask(:,:,1))/255;
figure;
imshow(abs(disparity).*Mask,[0 max_dis]);colorbar;

% Question 9:
% Comment the result obtained. Why do you think the disparity is not 
% well estimated in parts of the table?

