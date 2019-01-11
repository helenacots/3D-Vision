function [ disparity, cost, weight ] = stereo_computation( Il, Ir, min_dis, max_dis, wsize, cost_function, bilateral)
% The input parameters are 7:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
% - bilateral: flag that activates the bilateral weights or not

    [h,w,c] = size(Il);
    disparity=zeros(h,w);
    cost=zeros(h,w);
    
    stdC = 12;
    stdP = 17;
    w_step = floor(wsize/2);
    
    for i=1+w_step:h-w_step
        for j=1+w_step:w-w_step
            windowl = double(Il(i-w_step:i+w_step, j-w_step:j+w_step,:));
                        
            kmin = max(1+w_step, j+min_dis);
            kmax = min(w-w_step, j+max_dis);
                        
            [wh,ww,c] = size(windowl);
            weight = repmat(1/(wh*ww),wh,ww);
            
            % Bilateral Weights
            if(bilateral)
                 % TO COMPLETE
                
                 Np = sum(windowl,3);
                 P = zeros(wsize);
                 centerPoint = floor(wsize/2);
                 P(:,:) = Np(centerPoint,centerPoint);
                 C = (1/3)*abs(P-Np);
                 wcol = exp(-C/stdC);
                 
                 Q = [];
                 q = [];
                 
                 for n = i-w_step:i+w_step
                     for m = j-w_step:j+w_step
                         q = [n,m];
                         Q = [Q;q];
                     end
                 end
                 
                 G = sqrt(sum(([i,j]-Q).^2,2));
                 wpos = exp(-G/stdP);
                 wpos = (reshape(wpos,[wsize,wsize]))';
                 
                 weight = wcol.*wpos;
                 weight = weight/max(max(weight));

            end
            
            if strcmpi(cost_function,'SSD')                      
                error = Inf;
                for k = kmin:kmax
                   
                    windowr = double(Ir(i-w_step:i+w_step, k-w_step:k+w_step,:));
                    costSSD = sum(sum(sum(weight.*abs(windowl-windowr).^2)));
                    
                    if(costSSD<error)
                        error = costSSD;
                        kbest = k;
                    end
                end
                
            elseif strcmpi(cost_function,'NCC')  
                error = -Inf;
                for k = kmin:kmax
                    windowr = double(Ir(i-w_step:i+w_step, k-w_step:k+w_step,:));
                    term1 = windowl - sum(sum(sum(weight.*windowl)));
                    term2 = windowr - sum(sum(sum(weight.*windowr)));
                    sigmaIl = sqrt(sum(sum(sum(weight.*term1.^2))));
                    sigmaIr = sqrt(sum(sum(sum(weight.*term2.^2))));
                    costNCC = sum(sum(sum(weight.*term1.*term2)))/(sigmaIl*sigmaIr);
                    
                    if(costNCC>error)
                        error = costNCC;
                        kbest = k;
                    end
                end
                
            else
                error(strcat(cost_function,' not a valid matching cost name.'));
            end
            
            disparity(i, j) = kbest-j;
            cost(i, j) = error;
        end
    end
end      