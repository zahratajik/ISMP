function dp = Label_of_patch(edge_)
% Labeling of each patch.
%
% Inputs:
%   - edge_: Degree of node
%
% Outputs:
%   - dp: The matrix includes the total degree of patches, x, and y of the center of the patches.

edge_pos = edge_;
edge_pos(edge_pos < 0) = 0;

edge_neg = edge_;
edge_neg(edge_neg >= 0) = 0;

[L1, num1]=bwlabel(edge_pos,8);
for ij=1:num1 %ij counter of area image
    [r,c]=find(L1==ij); %length r refers to number of pixels of each magnetic patches
    s_pos = regionprops(L1,'centroid');
    cc_pos = s_pos(ij).Centroid;
    for y=1:length(r)
        degree_patch(y,1)=edge_(r(y),c(y));
        dpp(ij,2)= cc_pos(1,1); % x center
        dpp(ij,3)= cc_pos(1,2);% y center
        
    end
    dpp(ij,1)=sum(degree_patch(:,1)); %sum the degree f positive patch  
end



[L2, nunx]=bwlabel(edge_neg,8);
for ij=1:nunx %ij counter of area image
    [r,c]=find(L2==ij); %length r refers to number of pixels of each magnetic patches
    s_neg = regionprops(L2,'centroid');
    cc_neg = s_neg(ij).Centroid;
    for y=1:length(r)
        degree_patch(y,1)=edge_(r(y),c(y));
        dnp(ij,2)= cc_neg(1,1);
        dnp(ij,3)= cc_neg(1,2);
    end
    dnp(ij,1)=sum(degree_patch(:,1)); %degree of negitive patch
    
end

dp_ = vertcat(dpp, dnp);
dp = flip(sortrows(abs(dp_),1));
%dp= flip(dp);
end