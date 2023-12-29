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
<img src="https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a" alt="Identification" width="800">
</p>
<!-- ![Identification](https://github.com/zahratajik/ISMP/assets/75752814/018bc214-f4f2-4011-a155-7ece547ed62a) -->


