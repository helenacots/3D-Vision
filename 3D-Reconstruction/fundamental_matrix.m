% FUNDMATRIX - computes fundamental matrix from 8 or more points
%
% Function computes the fundamental matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.
%
% Usage:   F = fundmatrix(x1, x2)
%         
% Arguments:
%          x1, x2 - Two sets of corresponding 3xN set of homogeneous
%          points.
%         
% Returns:
%          F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
%

% Based on Peter Kovesi's code

function F = fundamental_matrix(x1, x2)

    npts = size(x1,2); 
    
    % Normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
    % normalise2dpts also ensures the scale parameter is 1.
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
    
    % Build the constraint matrix
    %Ai*f=0
    for i=1:npts
        
        A(i,:) = [x1(1,i)*x2(1,i), x2(1,i)*x1(2,i), x2(1,i), x2(2,i)*x1(1,i), x2(2,i)*x1(2,i), x2(2,i), x1(1,i), x1(2,i), 1];

    end
    
   
	[U,D,V] = svd(A,0); % use the economy decomposition
    %f is the las column of V --> V o V' ?? 
    f = V(:,size(V,2));
    
    %deberiaser0 = A*f; %molt proper a 0
    
    % Enforce constraint that fundamental matrix has rank 2 by performing
    % a svd and then reconstructing with the two largest singular values.
    F = [f(1) f(2) f(3);
         f(4) f(5) f(6);
         f(7) f(8) f(9)];
    [U,D,V] = svd(F); %% NO ENTENC PQ FAN AQUEST PAS
    D(3,3) = 0;    %bc we get the larger singular values s1 and s2
    F=U*D*V';
    
% % %     epipol1 = transpose(x2) * F;
% % %     epipol2 = F * x1; 
    
    %E = epipol2*epipol1
    % el rank ha de ser 2
    rankF = rank(F) ;
%     
%     % Denormalise
   F = T2'*F*T1;
    