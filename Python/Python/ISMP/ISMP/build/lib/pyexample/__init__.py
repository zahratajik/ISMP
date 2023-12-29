# Install Library
pip install networkx[all]

#Import Library
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat
import scipy.io
import networkx as nx


# Extraction of Connections

def strat_and_end_point_each_lines(im, n2, m2, i1, j1, i, j):

    """
    Calculating the Initial and Final Positions of Each Line:

    . im (array) = Image Magnetogram.
    . m2 and n2 (int) =  m2 and n2 are  Number of rows and columns, respectively.
    . i1 and j1 (array) = Indices of positive cells in magnetogram.
    . i and j (array) = Indices of negitive cells in magnetogram.

    Return :
    .  I (array) : Coordinates the beginning and Ending Positions of Lines.

    """
    
    numSamples = 2 * max(m2, n2) # number of points for creating a line
    im[w2] = np.abs(im[w2])
    I = np.zeros((len(i1), 2), dtype=object)
    I[0,:] = 0
    for k in range(len(i1)):
        p = 0
        end_ = [] #counter
        x1, y1 = j1[k], m2 - i1[k] + 1    #line start point coordinates
        for q in range(len(i)):
            x2, y2 = j[q], m2 - i[q] + 1    #line end point coordinates


            #Recognition of Pixels
            x = np.linspace(x1, x2, numSamples)
            y = np.linspace(y1, y2, numSamples)
            xy = np.round(np.column_stack([x,y]))
            dxy = np.abs(np.diff(xy, axis=0))
            duplicateRows = [0] + list(np.sum(dxy, axis=1) == 0)
            finalxy = xy[np.logical_not(duplicateRows)]


            # Coordinates of pixels laied on the line. Note, start and end points are excluded.
            finalx = finalxy[1:-1, 0]
            finaly = finalxy[1:-1, 1]


            # Indices of pixels laied on the line. Note, start and end points are excluded.
            ind_x1 = m2 - finaly + 1
            ind_x = ind_x1.astype(int)
            ind_y1 = finalx
            ind_y = ind_y1.astype(int)


            #Intensity of pixels laied on the line.
            intensity = np.zeros(len(ind_x))
            if len(ind_x) > 2:
               for ww in range(len(ind_x)):
                   intensity[ww]= im [ind_x[ww]-1, ind_y[ww]-1]


              # VG condition & adjacency matrix
               if (max(im[i1[k]-1, j1[k]-1], im[i[q]-1, j[q]-1]) > np.max(intensity)):
                  lin_ind_2 = np.ravel_multi_index ((i[q], j[q]), im.shape) # linear index of  end  point for adjacency matrix
                  end_.append(lin_ind_2)
                  p = p + 1
            else:
                lin_ind_2 = np.ravel_multi_index ((i[q], j[q]), im.shape) # linear index of  end  point for adjacency matrix
                end_.append(lin_ind_2)
                lin_ind_1 = np.ravel_multi_index((i1[k], j1[k]), im.shape)
        lin_ind_1 = np.ravel_multi_index((i1[k], j1[k]), im.shape) # linear index of start point for adjacency matrix

        if len(end_) != 0:
           I[k][0] = lin_ind_1
           I[k][1] = end_
    return I

# Weighted Graph's Adjacency Matrix
def create_adj_matrix(m2,n2, im, I, i1):

    """
    Calculating Adjacency Matrix :
    . m2 and n2 (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.
    . im (array) =  Image Magnetogram
    . I (array) = Coordinates the beginning and Ending Positions of Lines.
    . i1 (array)= Indices of Positive Cells
    Return:
    . adj (array)= Adjacency Matrix

    """

    N_pix = m2 * n2                    # number of pixel in each image
    adj = np.zeros((N_pix, N_pix)) # adjacancy matrix
    for k in range(len (sum(I[:,0] == '' for i in I))):
        sp = I[k][0] # start point linear index
        ep = np.asarray(I[k][1]) # end  point linear index
        i1, j1 = np.unravel_index(sp, im.shape) # start point subscripts
        adj[sp, ep] = im[i1, j1]
        for l in range(len(ep)):
            i2, j2 = np.unravel_index(ep[l], im.shape) # end point subscripts
            adj[ep[l], sp] = -im[i2, j2]
    return adj

# Degree of Nodes
def degree_of_node (adj, m2, n2):

    """
    Calculating Degree of Nodes:
    . adj (array) = Adjacency matrix
    . m2 and n2 (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.

    Return:
    . edge_ (array): Degree of Node.

    """
    edge = np.sum(adj, axis=1)
    edge_ = edge.reshape((m2, n2))
    return edge_

# Page-Rank
def page_rank (G, m2, n2):

    """
    Calculating Page-Rank:
    . G (Diagraph) = Graph
    . m2 and n2 (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.

    Return:
    . pr_ (array) = Page-Rank
    """
    pr=nx.pagerank(G, alpha=0.85)
    new_pr = list(pr.items())  #  convert page-rank to array
    pr_arr = np.array(new_pr)
    pr_ = pr_arr[:,1].reshape((m2, n2))
    return pr_

# Plot
def Plot (pr_,edge_,im1):

    """
    Plot :
    . pr_ (array) = Page-Rank
    . edge_ (array) = Degree of Node
    . im1 (array) = Image Magnetogram
    """
    fig, axs = plt.subplots(1, 3, figsize=(17,3))
    plt.gca().set_aspect('equal', adjustable='box')


    im1_cmap = plt.get_cmap('gray')
    im1_im = axs[0].imshow(im1, cmap=im1_cmap, vmin=-60, vmax=60)
    axs[0].set_xlabel('Solar-X', fontsize=18, fontname='Times New Roman')
    axs[0].set_ylabel('Solar-Y', fontsize=18, fontname='Times New Roman')
    axs[0].set_title('2013.11.01-01:19:13', fontsize=18, fontname='Times New Roman')
    axs[0].axis('on')
    axs[0].invert_yaxis()
    axs[0].tick_params(axis='both', which='major', labelsize=18)


    edge_cmap = ListedColormap(['black', 'gray', 'white'])
    edge_ = edge_ * 60
    edge_im = axs[1].imshow(edge_, cmap=edge_cmap, vmin=-60, vmax=60)
    axs[1].set_xlabel('Solar-X', fontsize=18, fontname='Times New Roman')
    axs[1].set_title('Map Degree of Node', fontsize=18, fontname='Times New Roman')
    axs[1].axis('on')
    axs[1].invert_yaxis()
    axs[1].tick_params(axis='both', which='major', labelsize=18)


    pr_cmap = plt.get_cmap('pink')
    pr_im = axs[2].imshow(pr_, cmap=pr_cmap)
    axs[2].set_xlabel('Solar-X', fontsize=18, fontname='Times New Roman')
    axs[2].set_title('Map of PageRank', fontsize=18, fontname='Times New Roman')
    axs[2].axis('on')
    axs[2].invert_yaxis()
    axs[2].tick_params(axis='both', which='major', labelsize=18)
    fig.colorbar(pr_im, ax=axs[2])

    plt.show()

    # Main
    im = fits.getdata("HMI20220117_235852_6173.fits") # Load data
    im1 = np.copy(im)
    m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively
    threshold = 25                     #threshold for magnetic field, in Gauss
    w1 = np.where(im > threshold)      #positive cells with values greater than a threshold
    w2 = np.where(im < -threshold)     #negative cells with values smaller than a threshold
    i1, j1 = w1[0], w1[1]              #indices of positive cells
    i, j = w2[0], w2[1]                #indices of negative cells
    I= strat_and_end_point_each_lines(im, n2, m2, i1, j1, i, j) # Extraction of connections
    adj = create_adj_matrix(m2,n2, im, I,i1)    # Weighted graph's adjacency matrix
    edge_= degree_of_node (adj, m2, n2)           # Degree of node
    G=nx.DiGraph(adj)                      # Graph
    pr_= page_rank (G, m2, n2)                  #Page-Rank
    Plot (pr_,edge_,im1)