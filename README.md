# ISMP Package (Identifying Solar Magnetic Patches)

The ISMP code package employs a complex network method to detect both large and small magnetic patches on the Sun's surface. The main routine and its corresponding functions are developed in MATLAB.
However, some parts of the algorithm are also available for Python users. In the future, the complete Python package will be available as well. 

## Description

The code includes the following functions:

**start_end_lines**: Calculating the initial and final positions of each line (Visibility Graph).

**adjacancy_matrix**: Calculating Adjacency Matrix.

**Label_of_patch**: Labeling of each patch.

**Rankdown**: Get the borders of the patches.

**pad_array**: Changing the size of the matrix used by the Rankdown function.

**Plot**: Plot Image.

**findDuplicatesAndMissing**: Find the number of duplicates in an array.

**Track**: Tracking magnetic patches.

## Examples

In the following example, we provide the identification and tracking magnetic patches for a given image. 
 
You can use these example images in the Examples Folder. 
Note: Run the follow code and test the examples with the codes provided below.

**Note**

For more information about this package, please refer to the following article.

[Tajik, Z., Farhang, N.,  Safari, H.,  & Wheatland, M. S., 2023](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=kIrSiKcAAAAJ&citation_for_view=kIrSiKcAAAAJ:2osOgNQ5qMEC)

## First-Identification
```matlab


clear all
close all
clc

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
Lp = Label_of_patch(edge_); % Label of each patch.
label_frist=(1:length(Lp)); % Labeling the magnetic patches of the first magnetogram
for i=1:length(label_frist)
    Lp(i).label=label_frist(i);
end
Plot(pr_,edge_,Image,Lp,pos_mask,neg_mask); 

```

<p align="center">
<img src="https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a" alt="Identification" width="1000">
</p>
<!-- ![Identification](https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a) -->

Second- Tracking

```matlab
directory = 'C:\Users\SSZ\Documents\GitHub\ISMP\Matlab\data-track';
final_matrix = Track(directory);
```

<p align="center">
<img src="https://github.com/zahratajik/ISMP/assets/75752814/89b430b7-d79b-42eb-95df-9b7ff1dd57a4" alt="Track" width="800">
</p>
<!-- ![Track](https://github.com/zahratajik/ISMP/assets/75752814/89b430b7-d79b-42eb-95df-9b7ff1dd57a4) -->





## Authors

[Zahra Tajik](https://scholar.google.com/citations?hl=en&user=kIrSiKcAAAAJ), [Nastaran Farhang](https://scholar.google.com/citations?hl=en&user=KGEB8dEAAAAJ) ,[Hossein Safari](https://scholar.google.com/citations?user=nCc1FV8AAAAJ&hl=en), [Michael S. Wheatland](https://scholar.google.com/citations?hl=en&user=m_oVye0AAAAJ)

