function   [img,mask]=rankdown(img,threshold)
img=pad_array(img,1);
if (threshold < 1.)
    threshold = .1*max(abs(img(:)));
end
nx=size(img,2);
ny=size(img,1);
fnx = fix(size(img,2));
mask=zeros(size(img,1),size(img,2));
offsets = [ 0, 1, ny+1, ny, ny-1, -1,  -ny-1, -ny, -ny+1];
nloop =2;
for q2 = 1:nloop
    if  q2 == 1
        tempimg = img;
    else
        tempimg = img*(-1);
    end
    
    % 1-D array of ranked pixel addresses
    %=====================================
    [sor, num] = sort(tempimg(:)); % highest-to-lowest
    ranks=flip(num);
    zer_find = find(tempimg < threshold);
    n_zer_find=length(zer_find);
    n_check = fnx*ny - n_zer_find; %only check above threshold
    
    tempmask = zeros(ny,nx);
    tempmask(ranks) = find(img(:)==img(:));
    
    if n_zer_find ~= 0
        tempmask(zer_find) = n_check;
    end
    
    current_label = 2;
    for q1 = 2:n_check
        
        % neighbor subarray
        %===================
        [maxnbr,maxpix] = max(tempimg(ranks(q1) + offsets ));
        nbrlabels = tempmask( ranks(q1) + offsets);
        
        % neighbor w/lowest label, or a new local maximum?
        %===================================================
        if length(maxpix)~= 0
            tempmask(ranks(q1)) = nbrlabels(maxpix);
        else
            tempmask(ranks(q1)) = current_label;
            current_label = current_label + 1;
        end
        
    end
    
    % offset 2nd iteration by max label from 1st
    %=============================================
    if  q2 == 2
        tempmask = tempmask + max(mask(:));
    end
    tempmask(zer_find) = 0 ;
    if (q2 == 2)
        tempmask = -tempmask;
    end
    mask = mask + tempmask;
    
end
img  =  img(2:end-1,2:end-1); 
mask = mask(2:end-1,2:end-1); %un-pad mask, too! same size = good
mask = fix(mask);

end
