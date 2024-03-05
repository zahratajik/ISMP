
clear all
close all
clc

%Identification

Image=fitsread('HMI20220117_235852_6173.fits'); %Load Magnetogram.
[nx,ny]=size(Image); % nx and ny are  Number of rows and columns, respectively.
threshold = 18; % threshold for magnetic field, in Gauss.
N_pix=nx*ny;  % Number of pixel in each image.
I =start_end_lines(Image, threshold); % Coordinates beginning and ending position of lines.
adj=adjacancy_matrix(Image,I); % Adjacency matrix.
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

%%

%Track

% Define the directory path
directory = 'C:\Users\SSZ\Documents\GitHub\ISMP\Matlab\data-track';
final_matrix = Track(directory); 