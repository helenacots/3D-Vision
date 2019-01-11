clear all;
close all;
addpath(genpath('.'));      



%% Test the triangulate function
% Use this code to validate that the function triangulate works properly
% Questions 1, 2 and 3

P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

% If triangulate.m is OK, this will be a matrix with very small numbers.
error = euclid(X_test) - euclid(X_trian);





%% Read template and images.
Irgb{1} = imread('images1/DSC_0495.JPG');
Irgb{2} = imread('images1/DSC_0503.JPG');
% you can start with this pair of images and then try with the rest of images available
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});


%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01); % threshold may be changed (see below)
    % start with a threshold value of 0.01, then decrease it to 0.005, and to 0.001.
    % this way you will get more correspondences and more 3D points
    % (sift function will be slowed down when you decrease the threshold)
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');


%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

vgg_gui_F(Irgb{1}, Irgb{2}, F');

%% Compute candidate camera matrices.

% Camera calibration matrix computed in practica3.m
K = [3044.19815016907 -4.01014339223785 1507.81572864665;
     0 3044.21428754626 1002.01971055401;
     0 0 1];

% Question 4
% Compute the Essential matrix from the Fundamental
% E = K2'*F*K1 , en aquest cas, com les dues vistes estan creades a partir
% de la mateixa càmera K2=K1.
E = K'*F*K;

[U,S,V] = svd(E);
% The SVD of E has several ambiguities.  In particular, U and V may be
% improper rotations in which case we need to change their sign.
if det(U) < 0
    U = -U;
end
if det(V) < 0
    V = -V;
end
W = [0 -1 0; 1 0 0; 0 0 1];
% E = SR 
%(S,R) = (U*Z*U',U*W*V')
% or (S,R) = (U*Z*U',U*W'*V')

% Question 5: Camera projection matrix P

% Assume that first camera is at the origin

P1 = K*[1 0 0 0; 0 1 0 0; 0 0 1 0];

% Question 6
% Rotation and translation for the 2nd camera can be computed from the SVD
% of E. There are 4 possible positions for the second camera.
Pc2 = {};
R1= U*W*V';
R2= U*W'*V';
u3 = U(:,size(U,2)); % last column of U;
Pc2{1} = K*[R1(:,:),u3(:)];
Pc2{2} = K*[R1(:,:),-u3(:)];
Pc2{3} = K*[R2(:,:),u3(:)];
Pc2{4} = K*[R2(:,:),-u3(:)];

figure;
plot_camera(P1,w,h);
plot_camera(Pc2{1},w,h);
plot_camera(Pc2{2},w,h);
plot_camera(Pc2{3},w,h);
plot_camera(Pc2{4},w,h);


%% Reconstruct structure
% Question 7: What does the following code do?
% --> choses the right solution from the 4 candidates
for i = 1:4
    X = triangulate(x1(:,1), x2(:,1), P1, Pc2{i}, [w h]);
    X = X / X(end);
    if P1(3,:) * X > 0 && Pc2{i}(3,:) * X > 0
        P2 = Pc2{i};
    end
end

% Question 8: Triangulate all matches
N = size(x1,2);
X = zeros(4,N);
for i = 1:N
    X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
end


%% Plot reconstructed points with colors
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
plot_camera(P1,w,h); % you may comment these two lines so as to better interact with the reconstructed 3D points
plot_camera(P2,w,h); % 
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end
axis equal;


%% Question 9: Compute reprojection error.
% The distance between the detected keypoints, x1 and x2, and the projection of the reconstructed
% points.

X(1:4,:) = X(1:4,:)./X(4,:);

% canvi de coordenades
%X1 = euclid(P1*X);
X1 = (P1*X);
% % X2 = euclid(P2*X);
X2 = (P2*X);
x1_h = homog(x1);
x2_h = homog(x2);

for i=1:N
reprojection_error(i) = (1/N)*sqrt((norm(P1*X(:,i)-x1_h(:,i))).^2 + (norm(P2*X(:,i) - x2_h(:,i))).^2);
end
figure;
hist(reprojection_error);




