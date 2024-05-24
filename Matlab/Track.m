function final_matrix = Track(directory)

%Tracking Magnetic Patches:
%Inputs:
%
%    . directory  =  data path.
%    

%Outputs:
%    . final_matrix = This matrix contains comprehensive information about the patches of each image, 
%such as the image name, patch labels, patch areas, center coordinates, edge matrix, 
%and false-positive index.


% Get a list of all .fits files in the directory
files = dir(fullfile(directory, '*.fits'));

n=1; % frist magnetogram
% Initialize variables
m1 = length(files);
fname = files(n,1).name;
Image1 = fitsread(fullfile(directory, fname));
threshold = 18; % threshold for magnetic field, in Gauss
[m2,n2] = size(Image1); % m2 and n2 are  Number of rows and columns, respectively
N_pix = m2*n2;  % number of pixel in each image

%Extraction adjacency and degree of node matrix

I = start_end_lines(Image1, threshold);
adj = adjacancy_matrix(Image1,I); % Adjacency matrix
edge = sum(adj,2); % Degree of node
edge_1 = reshape(edge,[m2, n2]); %Reshape degree of node to image size

Lp_1 = Label_of_patch(edge_1); % Label of Patches

label_frist=(1:length(Lp_1)); % Labeling the magnetic patches of the first magnetogram
for i=1:length(label_frist)
    Lp_1(i).label=label_frist(i);
end

[Image1,pos_mask1]=rankdown(edge_1,threshold); %Extract boundries
[Image1,neg_mask1]=rankdown(-edge_1,threshold); %Extract boundries

% imshow frist magnetogram with labels and contours
figure,subplot(1,2,1),imshow(edge_1,[-60 60])
hold on
axis on
axis xy
contour(pos_mask1,1,'LineWidth',2,'LineColor','r')
contour(neg_mask1,1,'LineWidth',2,'LineColor','r')
set(gca, 'FontName', 'Times New Romhisan', 'FontSize', 22)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22);
xlabel('Solar-X', 'FontSize', 22, 'FontName', 'Times New Roman')
ylabel('Solar-Y', 'FontSize', 22, 'FontName', 'Times New Roman')

title('Step 1')

for i=1:numel(Lp_1)
    text(round(Lp_1(i).xc), round(Lp_1(i).yc), sprintf('%d',Lp_1(i).label), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','color','b',...
        'FontSize',12,'FontName','Times New Roman',...
        'FontWeight','bold');
    hold on
end



final_matrix(n) = save_Info(fname, Lp_1);


for n=2:m1
    clear Image fname I new_label adj edge edge_2 Lp_2 common_rows PL repeated_1  repeated_2  new_label ind index
    jj=1;
    
    %Magnetogram
    fname = files(n,1).name;
    Image = fitsread(fullfile(directory,fname));
    
    %Extraction adjacency and degree of node matrix
    I =start_end_lines(Image, threshold);
    adj=adjacancy_matrix(Image,I); % Adjacency matrix
    edge = sum(adj,2); % Degree of node
    edge_2 =reshape(edge,[m2, n2]); %Reshape degree of node to image size
    Lp_2 = Label_of_patch(edge_2);% Label of patches
    
    if n==2
        
        % Examining the percentage of similarity and proximity of the centers of two patches
        for k = 1:sum(~cellfun(@isempty,{Lp_1.xc}))
            xyp_1(:,1) = Lp_1(k).x; % x all pixel in each patch for image 1
            xyp_1(:,2) = Lp_1(k).y; % y all pixel in each patch for image 1
            xyc_1(1,1)=Lp_1(k).xc; % x center
            xyc_1(1,2)=Lp_1(k).yc;% y center
            size_1=Lp_1(k).A; %Area
            for kk =  1:sum(~cellfun(@isempty,{Lp_2.xc}))
                xyp_2(:,1) = Lp_2(kk).x; % x all pixel in each patch for image 2
                xyp_2(:,2) = Lp_2(kk).y; % y all pixel in each patch for image 2
                xyc_2(1,1)=Lp_2(kk).xc; % x center
                xyc_2(1,2)=Lp_2(kk).yc; % y center
                size_2=Lp_2(kk).A; %Area
                D = sqrt(max(size_1,size_2));
                if  abs(xyc_1(1,1)-xyc_2(1,1)) < D
                    common_rows  = intersect(xyp_1(:, :), xyp_2(:, :), 'rows');
                    tresh = (length(common_rows(:,1))/length(xyp_1(:,1)))*100; % Similarity percentage of two patches
                    if tresh >= 30
                        Lp_2(kk).label=("O-"+""+num2str(Lp_1(k).label));
                        PL(jj,1) = k; %information of the patches' labels 
                        PL(jj,2) = kk; %information of the patches' labels 
                        jj=jj+1;
                    end
                end
                clear xyp_2 xyc_2 tresh size_2 r
            end
            clear xyp_1 xyc_1 size_1
        end
        
        %The emergence of the magnetic patch
        ind = find(strcmp({Lp_2.label}, 'n'));
        
        h=1;
         for p=1:sum(~cellfun(@isempty,{Lp_1.label}))
                old_label(h) = (Lp_1(p).label);
                h=h+1;
        end
        
        for i2=1:length(ind)
            newlabel_Array = {num2str(max(old_label)+i2)};
            Lp_2(ind(i2)).label = ("E-"+""+newlabel_Array);
        end
        
        
        %% Finding merge and fragmentation
         %% Finding merge and fragmentation
        repeated_1 = findDuplicates(PL(:,1)); % Find Fragment
        %Fragment
        if ~isempty(repeated_1)
            Lp_2 = Fragmentation(repeated_1, PL, Lp_2, old_label, ind,qs);
        end
        
     
        PL(any(ismember(PL,repeated_1), 2), :) = []; % remove the index corresponding with Fragment
        
       %Merge
        repeated_2 = findDuplicates(PL(:,2)); % Find Merge
        if repeated_2 ~= 0
          Lp_2 = Merge(repeated_2, PL, Lp_1, Lp_2, n);
        end
        
    else
        
        % Examining the percentage of similarity and proximity of the centers of two patches
        for k = 1:sum(~cellfun(@isempty,{Lp_1.xc}))
            xyp_1(:,1) = Lp_1(k).x; % x all pixel in each patch for image 1
            xyp_1(:,2) = Lp_1(k).y; % y all pixel in each patch for image 1
            xyc_1(1,1)=Lp_1(k).xc; % x center
            xyc_1(1,2)=Lp_1(k).yc;% y center
            size_1=Lp_1(k).A; %Area
            for kk =  1:sum(~cellfun(@isempty,{Lp_2.xc}))
                
                xyp_2(:,1) = Lp_2(kk).x; % x all pixel in each patch for image 2
                xyp_2(:,2) = Lp_2(kk).y; % y all pixel in each patch for image 2
                xyc_2(1,1)=Lp_2(kk).xc; % x center
                xyc_2(1,2)=Lp_2(kk).yc; % y center
                size_2=Lp_2(kk).A; %Area
                D = sqrt(max(size_1,size_2));
                
                if  abs(xyc_1(1,1)-xyc_2(1,1)) < D
                    common_rows  = intersect(xyp_1(:, :), xyp_2(:, :), 'rows');
                    tresh = (length(common_rows(:,1))/length(xyp_1(:,1)))*100; % Similarity percentage of two patches
                    if tresh >= 30
                        numbers = regexp (Lp_1(k).label, '[0-9]+', 'match');
                        if length(numbers(1,:))== 1
                            Lp_2(kk).label=("O-"+""+numbers);
                        else
                            Lp_2(kk).label=("O-"+""+numbers(1,2));
                        end
                        
                        PL(jj,1) = k;
                        PL(jj,2) = kk;
                        jj=jj+1;
                    end
                end
                clear xyp_2 xyc_2 tresh size_2 r
            end
            clear xyp_1 xyc_1 size_1
        end
        
        %The emergence of the magnetic patch
        ind = find(strcmp({Lp_2.label}, 'n'));
        
        h=1;
        for p=1:sum(~cellfun(@isempty,{Lp_1.label}))
            extrct_num=(regexp ((Lp_1(p).label), '[0-9]+', 'match'));
            if length(extrct_num)==1
                old_label(h) = str2double(string(extrct_num));
            else
                old_label(h) = str2double(string(extrct_num(1,2)));
            end
            h=h+1;
        end
        
        for i2=1:length(ind)
            newlabel_Array = {num2str(max(old_label)+i2)};
            Lp_2(ind(i2)).label = ("E-"+""+ newlabel_Array);%The emergence of the magnetic patch
        end
        
        
        %% Finding merge and fragmentation
        
         %Fragment
        repeated_1 = findDuplicates(PL(:,1)); % Find Fragment
        if ~isempty(repeated_1)
            Lp_2 = Fragmentation(repeated_1, PL, Lp_2, old_label, ind,qs);
        end
        
        PL(any(ismember(PL,repeated_1), 2), :) = [];
        % Merge
        repeated_2 = findDuplicates(PL(:,2)); % Find Merge
        if repeated_2 ~= 0
          Lp_2 = Merge(repeated_2, PL, Lp_1, Lp_2, n);
        end
        
    end
    clear yc xc A label
    
    % save information each magntogram
    final_matrix(n) = save_Info(fname, Lp_2);
    
    
    if n>=3 && n<m1-1
        % Find index of noise patches
        index = chek_noise_patch(final_matrix,n);
        if ~isempty(index)
            % Remove noise patches
            final_matrix = remove_noisie_patch(final_matrix,n, index);
        end
    end
 
    clear Lp_1
    Lp_1 = Lp_2;
end

% imshow latest magnetogram
edge_2=final_matrix(m1).edge;
label=final_matrix(m1).label;


[Image,pos_mask2]=rankdown(edge_2,threshold); %Extract boundries
[Image,neg_mask2]=rankdown(-edge_2,threshold); %Extract boundries
subplot(1,2,2), imshow(edge_2,[-60 60])
hold on
axis on
axis xy
set(gca, 'FontName', 'Times New Romhisan', 'FontSize', 22)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22);
xlabel('Solar-X', 'FontSize', 22, 'FontName', 'Times New Roman')
contour(pos_mask2,1,'LineWidth',2,'LineColor','r')
contour(neg_mask2,1,'LineWidth',2,'LineColor','r')
title(sprintf('Step %d', m1))
for i=1:numel(Lp_2)
    text(round(Lp_2(i).xc), round(Lp_2(i).yc), sprintf('%s',string(label(i))), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','color','b',...
        'FontSize',12,'FontName','Times New Roman',...
        'FontWeight','bold');
    hold on
end

end