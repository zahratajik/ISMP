# Install Library
#pip install networkx[all]

#Import Library
from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from scipy.io import savemat
from astropy.io import fits
import networkx as nx
import numpy as np
import scipy.io
import warnings

"""Extraction of Connections"""

def start_end_lines(im_, ny_, nx_,threshold):

    """
    Calculating the Initial and Final Positions of Each Line:
        
    Parameters
    ----------    

    . im (array) = Image Magnetogram.
    . nx and ny (int) =  m2 and n2 are  Number of rows and columns, respectively.
    . i1 and j1 (array) = Indices of positive cells in magnetogram.
    . i and j (array) = Indices of negitive cells in magnetogram.

    Return :
    .  I (array) : Coordinates the beginning and Ending Positions of Lines.

    """
    numSamples = 2 * max(nx_, ny_) # number of points for creating a line
    im_[w2_] = np.abs(im_[w2_])
    I = np.zeros((len(i1_), 2), dtype=object)
    I[0,:] = 0
    w1_ = np.where(im > threshold)      #positive cells with values greater than a threshold
    w2_ = np.where(im < -threshold)     #negative cells with values smaller than a threshold
    i1_, j1_ = w1_[0], w1_[1]              #indices of positive cells
    i_, j_ = w2_[0], w2_[1]                #indices of negative cells
    for k in range(len(i1_)):
        p = 0
        end_ = [] #counter
        x1, y1 = j1_[k], nx_ - i1_[k] + 1    #line start point coordinates
        for q in range(len(i_)):
            x2, y2 = j_[q], nx_ - i_[q] + 1    #line end point coordinates


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
            ind_x1 = nx_ - finaly + 1
            ind_x = ind_x1.astype(int)
            ind_y1 = finalx
            ind_y = ind_y1.astype(int)


            #Intensity of pixels laied on the line.
            intensity = np.zeros(len(ind_x))
            if len(ind_x) > 2:
               for ww in range(len(ind_x)):
                   intensity[ww]= im_ [ind_x[ww]-1, ind_y[ww]-1]


              # VG condition & adjacency matrix
               if (max(im_[i1_[k]-1, j1_[k]-1], im_[i_[q]-1, j_[q]-1]) > np.max(intensity)):
                  lin_ind_2 = np.ravel_multi_index ((i_[q], j_[q]), im_.shape) # linear index of  end  point for adjacency matrix
                  end_.append(lin_ind_2)
                  p = p + 1
            else:
                lin_ind_2 = np.ravel_multi_index ((i_[q], j_[q]), im_.shape) # linear index of  end  point for adjacency matrix
                end_.append(lin_ind_2)
                lin_ind_1 = np.ravel_multi_index((i1_[k], j1_[k]), im_.shape)
        lin_ind_1 = np.ravel_multi_index((i1_[k], j1_[k]), im_.shape) # linear index of start point for adjacency matrix

        if len(end_) != 0:
           I[k][0] = lin_ind_1
           I[k][1] = end_
    return I

"""Weighted Graph's Adjacency Matrix"""

def create_adj_matrix(m2_,n2_, im_, I_):

    """
    Calculating Adjacency Matrix :
    . m2 and n2 (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.
    . im (array) =  Image Magnetogram
    . I (array) = Coordinates the beginning and Final Positions of Lines

    Return:
    . adj (array)= Adjacency Matrix

    """

    N_pix = m2_ * n2_                    # number of pixel in each image
    adj_ = np.zeros((N_pix, N_pix)) # adjacancy matrix
    for k in range(len (sum(I_[:,0] == '' for i in I_))):
        sp = I_[k][0] # start point linear index
        ep = np.asarray(I_[k][1]) # end  point linear index
        i1_, j1_ = np.unravel_index(sp, im_.shape) # start point subscripts
        adj_[sp, ep] = im_[i1_, j1_]
        for l in range(len(ep)):
            i2, j2 = np.unravel_index(ep[l], im_.shape) # end point subscripts
            adj_[ep[l], sp] = -im_[i2, j2]
    return adj_

"""Degree of Nodes"""

def degree_of_node (adj_, nx_, ny_):

    """
    Calculating Degree of Nodes:
        
    Parameters
    ----------
    
    . adj (array) = Adjacency matrix
    . nx and ny (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.

    Return:
    . edge_ (array): Degree of Node.

    """
    edge = np.sum(adj_, axis=1) #Degree of node
    edge_ = edge.reshape((nx_, ny_)) #Reshape degree of node to image size.
    return edge_

"""Page-Rank"""

def page_rank (G_, nx_, ny_):

    """
    Calculating Page-Rank:
        
    Parameters
    ----------
    
    . G (Diagraph) = Graph
    . nx and ny (int) =  m2 and n2 are  Number of Rows and Columns, Respectively.

    Return:
    . pr_ (array) = Page-Rank
    """
    pr=nx.pagerank(G_, alpha=0.85)
    new_pr = list(pr.items())  #  convert page-rank to array
    pr_arr = np.array(new_pr)
    pr_ = pr_arr[:,1].reshape((nx_, ny_)) # Reshape Page-Rank to image size
    return pr_

""" Plot"""

def Plot (pr_,edge_,im1_):
    
    """
    Plot :
       
    Parameters
    ----------
    
    . pr_ (array) = Page-Rank
    . edge_ (array) = Degree of Node
    . im1 (array) = Image Magnetogram
    """
   
  
    fig, axs = plt.subplots(1, 3, figsize=(17,3))
    plt.gca().set_aspect('equal', adjustable='box')


    im1_cmap = plt.get_cmap('gray')
    im1_im = axs[0].imshow(im1_, cmap=im1_cmap, vmin=-60, vmax=60)
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


def run(FitsFile_ , threshold_=18):


    """
    
   
    Parameters
    ----------
    FitFile : Load Magnetogram.
     DESCRIPTION. The default load "HMI20220117_235852_6173.fits.".
    threshold (optional): Threshold for magnetic field, in Gauss.
        DESCRIPTION. The default is 18.


    """
im = fits.getdata("FitsFile.fits") # Load data
im1 = np.copy(im)
m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively
threshold = 18                     #threshold for magnetic field, in Gauss
I= start_end_lines(im, n2, m2, threshold) # Extraction of connections
adj=create_adj_matrix(m2,n2, im, I)   # Weighted graph's adjacency matrix
edge_= degree_of_node (adj, m2, n2)      # Degree of node
G=nx.DiGraph(adj)                        # Graph
pr_= page_rank (G, m2, n2)                  #Page-Rank
Plot (pr_,edge_,im1)