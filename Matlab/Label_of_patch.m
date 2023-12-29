function Lp = Label_of_patch(edge_)
% Calculating label of each patch.
%
% Inputs:
%   - edge_: Degree of node
%
% Outputs:
%   - dp: The matrix includes the total degree of each patch, x, and y of the center of the patches.

edge_pos = edge_;
edge_pos(edge_pos < 0) = 0;

edge_neg = edge_;
edge_neg(edge_neg >= 0) = 0;
%threshold downhill=25;
[L1, num1]=bwlabel(edge_pos,8);
for ij=1:num1 %ij counter of area image
    [y,x]=find(L1==ij); %length x refers to number of pixels of each magnetic patches
    
    s_pos = regionprops(L1,'centroid');
    cc_pos = s_pos(ij).Centroid;
    
    % Lpp = Label positive patches
    Lpp(ij).x= x; %x all pixel in each patch
    Lpp(ij).y= y; %y all pixel in each patch
    Lpp(ij).A= length(x); %Area
    Lpp(ij).xc= cc_pos(1,1); %x center
    Lpp(ij).yc= cc_pos(1,2);%y center
    
    
end

clear x y
[L2, num2]=bwlabel(edge_neg,8);
for ij=1:num2 %ij counter of area image
    [y,x]=find(L2==ij); %length r refers to number of pixels of each magnetic patches
    
    s_neg = regionprops(L2,'centroid');
    cc_neg = s_neg(ij).Centroid;
    
    % Label negative patches
    Lnp(ij).x= x; %x all pixel in each patch
    Lnp(ij).y= y; %y all pixel in each patch
    Lnp(ij).A= length(x);%Area
    Lnp(ij).xc= cc_neg(1,1); %x center
    Lnp(ij).yc= cc_neg(1,2);% y center
    
    
end


%The combination of two structures
Lp = [Lpp,Lnp];

%A selection of large patches
%pp=1;
%for ii=1:length(Lp1)
%   xyp = Lp1(ii).x;
 %   if length(xyp)>=5
%    Lp(pp).x= Lp1(ii).x; %x all pixel in each patch
 %   Lp(pp).y= Lp1(ii).y; %y all pixel in each patch
 %   Lp(pp).A= Lp1(ii).A;%Area
 %   Lp(pp).xc= Lp1(ii).xc; %x center
 %   Lp(pp).yc= Lp1(ii).yc;% y center 
  % pp=pp+1;
  %  end
%end

% Sort the structure according to values in descending order
[sortedColumn, sortIndex] = sort([Lp(:).xc],'ascend');

% Save the sorted output
Lp = Lp(sortIndex);

% label number
a=[1:length(Lp)];
for i=1:length(a)
    Lp(i).label= 'n';
end

end