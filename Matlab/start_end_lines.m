function I =start_end_lines(Image, threshold)

% Calculating the initial and final positions of each line:
%
%Inputs:
%    . Image  = Image Magnetogram.
%
%    . threshold = Threshold for magnetic field, in Gauss
%
%Outputs: :
%    .  I  : Coordinates beginning and ending position of lines.
[nx,ny]=size(Image); % nx and ny are  Number of rows and columns, respectively.
numSamples = 2 * max(nx, ny); % number of points for creating a line

w1 = find( Image > threshold );         % positive cells with values greater than a threshold
w2 = find( Image < -threshold );        % negative cells with values smaller than a threshold
[i1, j1] = ind2sub(size(Image), w1);    % indices of positive cells
[i, j] = ind2sub(size(Image), w2);      % indices of negative cells
Image(w2) = abs(Image(w2));
I(1).s = 0; I(1).e = 0;

for k = 1:length(i1)
    p = 1; end_ = []; % counter
    x1 = j1(k);    y1 = nx - i1(k) + 1;  % line start point coordinates
    for q = 1:length(i)
        x2 = j(q);  y2 = nx - i(q) + 1;  % line end point coordinates
        % Recognition of Pixels
        x = linspace(x1, x2, numSamples);
        y = linspace(y1, y2, numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0;sum(dxy, 2) == 0];
        %duplicateRows1 = [sum(dxy, 2) == 0];
        
        finalxy = xy(~duplicateRows,:);
        % Coordinates of pixels laied on the line. Note, start and end points are excluded.
        finalx = finalxy(2:end-1, 1);
        finaly = finalxy(2:end-1, 2);
        % Indices of pixels laied on the line. Note, start and end points are excluded.
        ind_x = nx - finaly + 1;
        ind_y = finalx;
        % Intensity of pixels laied on the line.
        for ww = 1:length(ind_x)
            intensity(ww) = Image(ind_x(ww), ind_y(ww));
        end
        % VG condition & adjacency matrix
        if ( max(Image(i1(k), j1(k)), Image(i(q), j(q))) > max(intensity) )
            lin_ind_2 = sub2ind(size(Image), i(q), j(q));   % linear index of  end  point for adjacency matrix
            end_(p) = lin_ind_2;
            p = p + 1;
        end
    end
    lin_ind_1 = sub2ind(size(Image), i1(k), j1(k)); % linear index of start point for adjacency matrix
    if  ~isempty(end_)
        I(k).s = lin_ind_1;
        I(k).e = end_;
    end
    
end

end



