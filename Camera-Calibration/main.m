clear all;
addpath(genpath('.'));


%% Read template and images.

%% Provided Images
T     = imreadnorm('images/template.JPG');
I{1}  = imreadnorm('images/DSC_0487.JPG');
I{2}  = imreadnorm('images/DSC_0488.JPG');
I{3}  = imreadnorm('images/DSC_0489.JPG');
I{4}  = imreadnorm('images/DSC_0491.JPG');
I{5}  = imreadnorm('images/DSC_0492.JPG');
N = length(I);

%% Our pattern images
% T     = imreadnorm('images1/template1.jpg');
% I{1}  = imreadnorm('images1/photo1.jpg');
% I{2}  = imreadnorm('images1/photo2.jpg');
% I{3}  = imreadnorm('images1/photo3.jpg');
% I{4}  = imreadnorm('images1/photo4.jpg');
% N = length(I);
 
%% Our Mafalda images
% T     = imreadnorm('images1/template2.jpg');
% I{1}  = imreadnorm('images1/maf1.jpg');
% I{2}  = imreadnorm('images1/maf2.jpg');
% I{3}  = imreadnorm('images1/maf3.jpg');
% I{4}  = imreadnorm('images1/maf4.jpg');
% N = length(I);

%% Optional: In case Matlab gets out of memory in your computer reduce the images size
% This will also accelerate computations

%I{1} = imresize(I{1}, 0.75);
%I{2} = imresize(I{2}, 0.75);
%I{3} = imresize(I{3}, 0.75);
% I{4} = imresize(I{4}, 0.75);
% I{5} = imresize(I{5}, 0.75);


%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(T, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end


%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x2 = points{i}(1:2, matches(2, :));
    H{i} = 0;
    [H{i}, inliers] = ransac_homography(homog(x1), homog(x2), 3, 1000);

    % Plot inliers.
    figure;
    plotmatches(T, I{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
    %vgg_gui_H(T, I{i}, H{i});
end


%% Compute the Image of the Absolute Conic

% Build the linear system of equations.
A = zeros(2 * N, 6);
for i = 1:N
    h = H{i};
    A(2*i - 1,:) = [h(1,1) * h(1,2), ...
                    h(1,1) * h(2,2) + h(2,1) * h(1,2), ...
                    h(1,1) * h(3,2) + h(3,1) * h(1,2), ...
                    h(2,1) * h(2,2), ...
                    h(2,1) * h(3,2) + h(3,1) * h(2,2), ...
                    h(3,1) * h(3,2)];
    A(2*i,:) = [h(1,1) * h(1,1), ...
                h(1,1) * h(2,1) + h(2,1) * h(1,1), ...
                h(1,1) * h(3,1) + h(3,1) * h(1,1) , ...
                h(2,1) * h(2,1), ...
                h(2,1) * h(3,1) + h(3,1) * h(2,1), ...
                h(3,1) * h(3,1)] ...
                -  [h(1,2) * h(1,2), ...
                h(1,2) * h(2,2) + h(2,2) * h(1,2), ...
                h(1,2) * h(3,2) + h(3,2) * h(1,2), ...
                h(2,2) * h(2,2), ...
                h(2,2) * h(3,2) + h(3,2) * h(2,2), ...
                h(3,2) * h(3,2)];
end

% Question 1
% Write the code to compute the IAC w by solving the above linear system. 

[U, S, V] = svd(A);
omega = V(:,length(V));
w = [omega(1),omega(2),omega(3);
            omega(2),omega(4),omega(5);
            omega(3),omega(5),omega(6)];
        
%% Recover the camera calibration.
% Question 2
% Compute the camera calibration matrix, K, from the IAC, w, using the
% Cholesky descomposition.

%K contains the internal parameters of the camera

K = inv(chol(w));
Knorm = K/K(3,3);

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure; hold;
for i = 1:N
    % Question 7
    
    h = H{i};
    r1 = inv(K)*h(:,1);         %#ok<MINV>
    r2 = inv(K)*h(:,2);           %#ok<MINV>
    t{i} = inv(K)*h(:,3);                 %#ok<MINV>

    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 3008, 2000, 200);
end

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Question 9: Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 3008, 2000, 200);
for i = 1:N
    
    vgg_scatter_plot([R{i}*p1+t{i}, R{i}*p2+t{i}, R{i}*p3+t{i}, R{i}*p4+t{i}, R{i}*p1+t{i}], 'r');
    
end


%% Question 11: Plot some 3D points on every camera.
[Th, Tw] = size(T);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';
% X are the coordinates of the vertices of a cube placed in the middle of
% the calibration target.
X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));
X = homog(X); %homogenous coordinates

for i = 1:N
    figure; colormap(gray);
    imagesc(I{i});
    hold on;
    x = P{i}*X;
    x=euclid(x);
    vgg_scatter_plot(x, 'c');
end

%% Question 12: Repeat all the process using the 5 provided images

%% Question 13: Repat all the process with your own images
    % IT WORKS!!
