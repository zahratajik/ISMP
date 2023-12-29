
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
dp = Label_of_patch(edge_); % Label of each patch.
Plot(pr_,edge_,Image,dp,pos_mask,neg_mask); 