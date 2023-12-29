function [index]=chek_noise_patch(final_matrix,n)

Middle_label = final_matrix(n-1).label;
Previous_label= final_matrix(n-2).label;
Next_label= final_matrix(n).label;
Middle_Area = final_matrix(n-1).Area;

%label middle image
Mid_Numbers= cellfun(@(x) str2double(regexp(x, '\d+', 'match')), Middle_label, 'UniformOutput', false);

for i = 1:length(Mid_Numbers)
    if length(Mid_Numbers{i}) == 1 %Remove small fragmented patches
        Mid_Numbers1(i) = Mid_Numbers{i};
    else
        Middle_Area(i)=0;
         
%     elseif length(Mid_Numbers{i}) == 2
%         Mid_Numbers1(i) = Mid_Numbers{i}(2);
    end
    clear w
end
w= find(Middle_Area==0);
Middle_Area(w1)=[];

% label next image
Next_Numbers= cellfun(@(x) str2double(regexp(x, '\d+', 'match')), Next_label, 'UniformOutput', false);

for i = 1:length(Next_Numbers)
    if length(Next_Numbers{i}) == 1 %Remove small fragmented patches
        Next_Numbers1(i) =Next_Numbers{i};
%     elseif length(Next_Numbers{i}) == 2
%         Next_Numbers1(i) = Next_Numbers{i}(2);
    end
end

if n-2>1
    
    % label pervious image
    prev_Numbers= cellfun(@(x) str2double(regexp(x, '\d+', 'match')), Previous_label, 'UniformOutput', false);
    
    for i = 1:length(prev_Numbers)
        if length(prev_Numbers{i}) == 1 %Remove small fragmented patches 
            prev_Numbers1(i) =prev_Numbers{i};
%         elseif length(prev_Numbers{i}) == 2
%             prev_Numbers1(i) = prev_Numbers{i}(2);
        end
    end
    
else
    prev_Numbers1 = Previous_label; 
end

w_small_area = find(Area < 3); % find patches with area 1 or 2 pixels
y=1;
for kk=1:length(w_small_area)
    w1 = find( prev_Numbers1 == Mid_Numbers1(w_small_area(kk)) ); % Searching  the label of the small patches of the middle image in the previous image
    if isempty(w1)
        w2 = find( Next_Numbers1 ==  Mid_Numbers1(w_small_area(kk)) ); %Searching  the label of the small patches of the middle image in the previous image
        if isempty(w2) 
          index(y)=w_small_area(kk); % Saving index of small patches that appeared only in the middle frame. 
          y=y+1;
        end
    end
end
end
