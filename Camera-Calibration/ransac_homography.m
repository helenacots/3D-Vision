function [H, idx_inliers] = ransac_homography(x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);

% ransac
it = 0;
best_inliers = [];
while it < max_it
    
    points = randomsample(Npoints, 4); 
    H = homography2d(x1(:,points), x2(:,points));
    inliers = compute_inliers(H, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    it = it + 1;

end

% compute H from all the inliers
H = homography2d(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(H, x1, x2, th)
    % Check that H is invertible
    if abs(log(cond(H))) > 15
        idx_inliers = [];
        return
    end
    
    % transformed points (in both directions)
    Hx1 = H * x1;
    Hix2 = inv(H) * x2;
    
    % normalise homogeneous coordinates (third coordinate to 1)
    x1 = normalise(x1);
    x2 = normalise(x2);     
    Hx1 = normalise(Hx1);
    Hix2 = normalise(Hix2); 
    
    % compute the symmetric geometric error
    d2 = sum((x1 - Hix2).^2) + sum((x2 - Hx1).^2);
    idx_inliers = find(d2 < th.^2);


function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat
