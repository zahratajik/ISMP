# ISMP Package (Identifying Solar Magnetic Patches)

The ISMP code package employs a complex network method to detect both large and small magnetic 
patches on the Sun's surface. The main routine and its corresponding functions are developed in MATLAB.
However, some parts of the algorithm are also available for Python users. In the future, the complete Python 
the package will be available as well. 

## Description

The code includes the following functions:

**start_end_lines**: Calculating the initial and final positions of each line (Visibility Graph).

**adjacancy_matrix**: Calculating Adjacency Matrix.

**Label_of_patch**: Labeling of each patch.

**Rankdown**: Get the borders of the patches.

**pad_array**: Changing the size of the matrix used by the Rankdown function.

**Plot**: Plot Image.

## Examples

In the following example, we provide the identification and tracking magnetic patches for a given image. 
 
You can use these example images in the Examples Folder. 
Note: Run the follow code and test the examples with the codes provided below.



## First-Identification
```matlab

clear all
close all
clc
Image=image;
Image=fitsread('HMI20220117_235852_6173.fits'); %Load Magnetogram.
threshold = 18; % threshold for magnetic field, in Gauss.
[nx,ny]=size(Image); % nx and ny are  Number of rows and columns, respectively.
N_pix=nx*ny;  % Number of pixel in each image.
I =start_end_lines(nx, ny, Image, threshold); % Coordinates beginning and ending position of lines.
adj=adjacancy_matrix(Image,I,nx,ny); % Adjacency matrix.
edge = sum(adj,2); % Degree of node.
edge_ =reshape(edge,[nx, ny]); %Reshape degree of node to image size.
G = graph(adj~=0); %Graph.
pr = centrality(G,'pagerank','FollowProbability',0.85); % Page-Rank.
pr_ =reshape(pr,[nx, ny]); %Reshape Page Rank to image size.
[img,pos_mask]=rankdown(edge_,threshold); %Extract boundries.
[img,neg_mask]=rankdown(-edge_,threshold); %Extract boundries.
dp = Label_of_patch(edge_); % Label of each patch.
Plot(pr_,edge_,Image,dp,pos_mask,neg_mask); 
```

<p align="center">
<img src="https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a" alt="Identification" width="1000">
</p>
<!-- ![Identification](https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a) -->

Second- Tracking

```matlab
clear all
clc
close all

% Define the directory path
directory = 'E:\file_PHD\HMI\crop_data\ARF2022_2990\extract';

% Get a list of all .fits files in the directory
files = dir(fullfile(directory, '*.fits'));

% Initialize variables
m1 = length(files);
fname = files(1,1).name;
Image1 = fitsread(fullfile(directory, fname));
threshold = 18; % threshold for magnetic field, in Gauss
[m2,n2] = size(Image1); % m2 and n2 are  Number of rows and columns, respectively
N_pix = m2*n2;  % number of pixel in each image

%Extraction adjacency and degree of node matrix

I = start_end_lines(m2, n2, Image1, threshold);
adj = adjacancy_matrix(Image1,I,m2,n2); % Adjacency matrix
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
figure(10),subplot(1,2,1),imshow(edge_1,[-60 60])
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

n=1; % frist magnetogram
for i=1:length(Lp_1)
    label(i)=Lp_1(i).label;
    A(i)= Lp_1(i).A;
    xc(i)=Lp_1(i).xc;
    yc(i)=Lp_1(i).yc;
end

final_matrix(n).imname = fname;
final_matrix(n).label = label;
final_matrix(n).Area = A;
final_matrix(n).xc = xc;
final_matrix(n).yc = yc;
final_matrix(n).edge=edge_1;



for n=2:10
    clear Image fname I new_label adj edge edge_2 Lp_2 common_rows PL repeated_1  repeated_2  new_label ind index
    jj=1;
    
    %Magnetogram
    fname = files(n,1).name;
    Image = fitsread(fullfile(directory,fname));
    
    %Extraction adjacency and degree of node matrix
    I =start_end_lines(m2, n2, Image, threshold);
    adj=adjacancy_matrix(Image,I,m2,n2); % Adjacency matrix
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
        
        % Finding new patches and labeling
        ind = find(cellfun(@(x) strcmp(x, 'n'), {Lp_2.label})); 
        
        h=1;
         for p=1:sum(~cellfun(@isempty,{Lp_1.label}))
                new_label(h) = (Lp_1(p).label);
                h=h+1;
        end
        
        for i2=1:length(ind)
            newlabel_Array = {num2str(max(new_label)+i2)};
            Lp_2(ind(i2)).label = ("E-"+""+newlabel_Array);
        end
        
        
        %% Finding merge and fragmentation
        repeated_1 = findDuplicatesAndMissing(PL(:,1)); % Find Fragment
        repeated_2 = findDuplicatesAndMissing(PL(:,2)); % Find Merge
        
        %Fragmentation
        if repeated_1 ~= 0
            for t=1:length(repeated_1)
                w1 (:,1) = find( PL(:,1) == repeated_1(t));
                fragment(:,1) = PL(w1,2);
                for j1=1:length(fragment(:,1))
                    Lp_2(fragment(j1,1)).label = sprintf("F%d-%d",j1,(max(new_label)+length(ind))+ j1);
                end
                
                clear fragment w1
            end
        end
        
        %Merge
        clear w1 Area w_max new_label2 w_A ww
        if repeated_2 ~= 0
            for t=1:length(repeated_2)
                ww(:,1) = find( PL(:,2) == repeated_2(t)); % Find duplicate numbers in step new
                merge(:,1) = PL(ww,1); % equivalent to the counters in step old
                for j2=1:length(merge(:,1))
                    Area(j2,1)= Lp_1(merge(j2,1)).A; % find Area of each patch
                end
                w_max = find( Area(:,1) == max(Area(:,1))); % maximum Area
                Lp_2(repeated_2(t)).label = ("M-"+""+ num2str(Lp_1(merge(w_max(1))).label));
                clear w_max ww merge Area
            end
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
        
        % Finding new patches and labeling
        ind = find(cellfun(@(x) strcmp(x, 'n'), {Lp_2.label}));
        
        h=1;
        for p=1:sum(~cellfun(@isempty,{Lp_1.label}))
                aaa=(regexp ((Lp_1(p).label), '[0-9]+', 'match'));
                if length(aaa)==1
                new_label(h) = str2num(string(aaa));
                else 
                new_label(h) = str2num(string(aaa(1,2))); 
                end
                h=h+1;
        end
        
        for i2=1:length(ind)
            newlabel_Array = {num2str(max(new_label)+i2)};
            Lp_2(ind(i2)).label = ("E-"+""+newlabel_Array);
        end
        
        
        %% Finding merge and fragmentation
        
        repeated_1 = findDuplicatesAndMissing(PL(:,1)); % find fragment
        repeated_2 = findDuplicatesAndMissing(PL(:,2)); % find collision
        
        %Fragmentation
        
        if repeated_1 ~= 0
            for t=1:length(repeated_1)
                w1 (:,1) = find( PL(:,1) == repeated_1(t));
                fragment(:,1) = PL(w1,2);
                for j1=1:length(fragment(:,1))
                    Lp_2(fragment(j1,1)).label = sprintf("F%d-%d",j1,(max(new_label)+length(ind))+ j1);
                    
                    
                end
                clear fragment w1 numbers
            end
        end
        clear fragment
        % Merge
        
        clear w1 Area w_max w_A ww
        
        if repeated_2 ~= 0
            for t=1:length(repeated_2)
                ww(:,1) = find( PL(:,2) == repeated_2(t)); % Find duplicate numbers in step new
                merge(:,1) = PL(ww,1); % equivalent to the counters in step old
                for j2=1:length(merge(:,1))
                    Area(j2,1)= Lp_1(merge(j2,1)).A; % find Area of each patch
                end
                w_max = find( Area(:,1) == max(Area(:,1))); % maximum Area
                
                numbers = regexp ((Lp_1(merge(w_max(1))).label), '[0-9]+', 'match');
                
                if length(numbers(1,:))==1
                    Lp_2(repeated_2(t)).label = ("M-" + "" + numbers);
                else
                    
                    Lp_2(repeated_2(t)).label = ("M-"+""+numbers(1,2));
                end
                clear w_max ww merge Area
            end
            
        end
        
    end
    clear yc xc A label
    
    % save information each magntogram
    label = cell(length(Lp_2), 1);
    for i=1:length(Lp_2)
        label{i} = Lp_2(i).label;
        A(i)=Lp_2(i).A;
        xc(i)=Lp_2(i).xc;
        yc(i)=Lp_2(i).yc;
    end
    
    final_matrix(n).imname = fname;
    final_matrix(n).label = label;
    final_matrix(n).Area = A;
    final_matrix(n).xc = xc;
    final_matrix(n).yc = yc;
    final_matrix(n).edge=edge_2;

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
edge_2=final_matrix(10).edge;
label=final_matrix(10).label;

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
```

<p align="center">
<img src="https://github.com/zahratajik/ISMP/assets/75752814/89b430b7-d79b-42eb-95df-9b7ff1dd57a4" alt="Track" width="800">
</p>
<!-- ![Track](https://github.com/zahratajik/ISMP/assets/75752814/89b430b7-d79b-42eb-95df-9b7ff1dd57a4) -->





## Authors

[Zahra Tajik](https://scholar.google.com/citations?hl=en&user=kIrSiKcAAAAJ), [Nastaran Farhang](https://scholar.google.com/citations?hl=en&user=KGEB8dEAAAAJ) , [Michael S. Wheatland](https://scholar.google.com/citations?hl=en&user=m_oVye0AAAAJ)

