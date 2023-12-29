function adj = adjacancy_matrix(Image,I,nx,ny)

%Calculating Adjacency Matrix :
%Inputs:
%    . nx and ny  =  nx and ny are  Number of Rows and Columns, Respectively.
%    . Image  =  Imageage Magnetogram
%    . I  = Coordinates the beginning and ending positions of lines.

%Outputs:
%    . adj = Adjacency Matrix

N_pix=nx*ny; % number of pixel in each Imageage
adj = zeros(N_pix,  N_pix); % adjacancy matrix
for k = 1:sum(~cellfun(@isempty,{I.s}))
    sp = I(k).s; % start point linear index
    ep = I(k).e; %  end  point linear index
    [i1, j1] = ind2sub(size(Image),sp); % start point subscripts
    adj (sp, ep(:)) = Image (i1, j1);
    for l = 1:length(ep)
        [i2, j2] = ind2sub(size(Image),ep(l)); % end point subscripts
        adj (ep(l), sp) = Image (i2, j2);
    end
end

