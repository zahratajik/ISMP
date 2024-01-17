# Install Library
# pip install networkx[all]
# pip install scikit-image

"""Import Library"""

from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import ListedColormap
from skimage.measure import regionprops
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import ndimage, io
from astropy.io import fits
import networkx as nx
import numpy as np
import scipy.io
import warnings



"""Extraction of Connections"""

def start_end_lines(im_, threshold_):

    """
    Calculating the Initial and Final Positions of Each Line:

    . im (array) = Image Magnetogram.
    . threshold_ (int) =  threshold_.
    
    output :
    .  I (array) : Coordinates the beginning and Ending Positions of Lines.

    """
    m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively
    w1_ = np.where(im_ > threshold_)      #positive cells with values greater than a threshold
    w2_ = np.where(im_ < -threshold_)     #negative cells with values smaller than a threshold
    i1_, j1_ = w1_[0], w1_[1]              #indices of positive cells
    i_, j_ = w2_[0], w2_[1]                #indices of negative cells
    numSamples = 2 * max(m2, n2) # number of points for creating a line
    im_[w2_] = np.abs(im_[w2_])
    I = np.zeros((len(i1_), 2), dtype=object)
    I[0,:] = 0
    for k in range(len(i1_)):
        p = 0
        end_ = [] #counter
        x1, y1 = j1_[k], m2 - i1_[k] + 1    #line start point coordinates
        for q in range(len(i_)):
            x2, y2 = j_[q], m2 - i_[q] + 1    #line end point coordinates


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


def create_adj_matrix( im_, I_):
    """
    Calculating Adjacency Matrix :

    . im (array) =  Image Magnetogram
    . I (array) = Coordinates the beginning and Final Positions of Lines

    output:
    . adj (array)= Adjacency Matrix

    """
    m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively.
    N_pix = m2 * n2                    # number of pixel in each image
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

def degree_of_node(adj_,im_):
    """
    Calculating Degree of Nodes:
    . adj (array) = Adjacency matrix.
    . im (array) =  Image Magnetogram
    output:
    . edge_ (array): Degree of Node.

    """
    m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively
    edge = np.sum(adj_, axis=1) #Degree of node
    edge_ = edge.reshape((m2, n2)) #Reshape degree of node to image size.
    return edge_


"""Page-Rank"""

def page_rank(G_,im_):
    """
    Calculating Page-Rank:
    . G (Diagraph) = Graph
    . im (array) =  Image Magnetogram.

    Return:
    . pr_ (array) = Page-Rank
    """
    m2, n2 = im.shape                  # m2 and n2 are  number of rows and columns, respectively
    pr=nx.pagerank(G_, alpha=0.85)
    new_pr = list(pr.items())  #  convert page-rank to array
    pr_arr = np.array(new_pr)
    pr_ = pr_arr[:,1].reshape((m2, n2)) # Reshape Page-Rank to image size
    return pr_


"""Label of patches"""

def label_of_patch(edge):

    """
      Determining various information for each patch, such as the center coordinate, area,
      coordinates of pixels within the patch, and labeling:
      input:

      . edge (array) = Degree of Node.
      output:

      . Lp (list) =  Includes  center coordinate, area, coordinates of pixels within the patch, and labeling.
    """
    edge_pos = np.maximum(edge, 0)
    edge_neg = np.minimum(edge, 0)
    threshold_value = 0  # Adjust this value based on your image
    binary_image_pos = (edge_pos > threshold_value).astype(int)

    # Perform connected component labeling
    L1, num1 = ndimage.label(binary_image_pos,structure=np.ones((3,3)))

    # Initialize a list to store properties of each patch
    LPP = []

    # Iterate through each labeled patch
    for ij in range(1, num1 + 1):
        # Find the coordinates of pixels in the current patch
        y, x = np.where(L1 == ij)

        # Calculate centroid
        s_pos = regionprops(L1)
        cc_pos = s_pos[ij - 1].centroid

        # Store properties in the patches_properties list
        pos_patch = {
            'label': ij,
            'x': x.tolist(),  # Convert to list for better compatibility
            'y': y.tolist(),
            'area': len(x),
            'center_x': cc_pos[1],
            'center_y': cc_pos[0]
        }

        LPP.append(pos_patch)

    binary_image_neg = (np.abs(edge_neg)> threshold_value).astype(int)

    # Perform connected component labeling
    L2, num2 = ndimage.label(binary_image_neg,structure=np.ones((3,3)))

    # Initialize a list to store properties of each patch
    LNP = []

    # Iterate through each labeled patch
    for ij in range(1, num2 + 1):
        # Find the coordinates of pixels in the current patch
        y, x = np.where(L2 == ij)

        # Calculate centroid
        s_neg = regionprops(L2)
        cc_neg = s_neg[ij - 1].centroid

        # Store properties in the patches_properties list
        neg_patch = {
            'label': ij,
            'x': x.tolist(),  # Convert to list for better compatibility
            'y': y.tolist(),
            'area': len(x),
            'center_x': cc_neg[1],
            'center_y': cc_neg[0]
        }
        LNP.append(neg_patch)

    LP=LPP+LNP


    # Sort the structure according to values in ascending order
    LP.sort(key=lambda patch: patch['center_x'])

    # Label number
    for i, patch in enumerate(LP):
      patch['label'] = i+1
    return LP

""" pad_array"""

def pad_array(arr, all):

   """
      Changing the size of the matrix:

      input:
      . arr (array) = Degree of Node.
      . all (int)= defult 1.

      output:
      . arr (array) = Matrix with new size.
    """
   if all == 1:
      arr1 = np.zeros((arr.shape[0]+1, arr.shape[1]))  #adds a row to top
      arr1[1:, :] = arr
      arr = arr1
   if all == 1:
      arr1 = np.zeros((arr.shape[0]+1, arr.shape[1])) #adds a row to bottom
      arr1[:-1, :] = arr
      arr = arr1
   if all == 1:
      arr1 = np.zeros((arr.shape[0], arr.shape[1]+1)) #adds col before
      arr1[:, 1:] = arr
      arr = arr1
   if all == 1:
      arr1 = np.zeros((arr.shape[0], arr.shape[1]+1)) # adds col after
      arr1[:, :-1] = arr
      arr = arr1
   return(arr)


"""rankdown"""

def rankdown(Image, threshold):

  """
      Get the borders of the patches:

      input:
      . Image (array) = Magnrtogram.
      . threshold (int)= Threshold.

      output:
      . mask (array): Boundaries matrix.

  """
  img = pad_array (Image,1)
  if threshold < 1:
    threshold = 0.1 * max(abs(img.ravel()))
  nx = img.shape[1]
  ny = img.shape[0]
  fnx = np.fix (img.shape[1])
  mask = np.zeros((img.shape[0], img.shape[1]))
  offsets = [0, 1, ny+1, ny, ny-1, -1, -ny-1, -ny, -ny+1]
  nloop = 1
  for q2 in range(0, nloop):
    tempimg= np.zeros((img.shape[0], img.shape[1]))
    if q2 == 0:
        tempimg = img.copy ()
    else:
        tempimg = img * (-1)
    sor, num = np.sort(tempimg.flatten()), np.argsort(tempimg.flatten())
    ranks = np.flip(num)
    zer_find = np.where(tempimg < threshold)
    n_zer_find = len (zer_find[0])
    n_check = int(fnx * ny - n_zer_find)
    tempmask = np.zeros((ny, nx))
    ranks= ranks.astype(int)
    tempmask.ravel()[ranks] = np.where(img.ravel() == img.ravel())
    if n_zer_find != 0:
      tempmask[zer_find[0],zer_find[1]] = n_check
    current_label = 1
    for q1 in range(1, n_check):
      q1=1
      maxnbr, maxpix = np.max(tempimg.ravel()[ranks[q1] + offsets]), np.argmax(tempimg.ravel()[ranks[q1] + offsets])
      nbrlabels = tempmask.ravel()[ranks[q1] + offsets]
      if len(str(maxpix)) != 0:
        tempmask.ravel()[ranks[q1]] = nbrlabels[maxpix]
      else:
        tempmask.ravel()[ranks[q1]] = current_label
        current_label += 1
      if q2 == 1:
        tempmask = tempmask + np.max(mask)
      tempmask[zer_find] = 0
      if q2 == 1:
          tempmask = -tempmask
      mask = mask + tempmask
  img = img[1:-1, 1:-1]
  mask = mask[1:-1, 1:-1]
  mask = np.fix(mask)
  return mask, img


""" Plot"""

def Plot (pr_,edge_,im1,LP):

    """
    Plot :
    input:

    . pr_ (array) = Page-Rank.
    . edge_ (array) = Degree of Node.
    . im1 (array) = Image Magnetogram.

    output:
    . Image
    """
    fig, axs = plt.subplots(1, 3, figsize=(12,3))
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
    axs[1].contour(mask1, levels=[1], linewidths=2, colors='r')
    axs[1].contour(mask2, levels=[1], linewidths=2, colors='r')
    # Assuming patches_properties is a list of dictionaries
    b1 = []
    c1 = []
    for i in LP:
        b1.append(i['center_x'])
        c1.append(i['center_y'])
        a1=(i['label'])
    b1 = np.array(b1)
    c1 = np.array(c1)
    for i in range(0,len(LP)):
        axs[1].text(round(b1[i]), round(c1[i]),str(i+1),
             horizontalalignment='center', verticalalignment='center',
             color='b', fontsize=8, fontname='Times New Roman', fontweight='bold')


    pr_cmap = plt.get_cmap('pink')
    pr_im = axs[2].imshow(pr_, cmap=pr_cmap)
    axs[2].set_xlabel('Solar-X', fontsize=18, fontname='Times New Roman')
    axs[2].set_title('Map of PageRank', fontsize=18, fontname='Times New Roman')
    axs[2].axis('on')
    axs[2].invert_yaxis()
    axs[2].tick_params(axis='both', which='major', labelsize=18)
    fig.colorbar(pr_im, ax=axs[2])

    plt.show()


def run(FitsFile, threshold_=18):

    """
    Perform various operations on a FITS file.
    
    Parameters
    ----------
    FitsFile : str, optional
    threshold_ : float, optional
        Threshold for magnetic field, in Gauss. The default is 18.
    """
    
    im = fits.getdata(FitsFile)  # Load data.
    im_copy = np.copy(im)
    m2, n2 = im.shape  # m2 and n2 are the number of rows and columns, respectively.
    threshold_ = 18
    I = start_end_lines(im, threshold_)  # Extraction of connections.
    adj = create_adj_matrix(m2, n2, im, I)  # Weighted graph's adjacency matrix.
    edge_ = degree_of_node(adj, m2, n2)  # Degree of node.
    G = nx.DiGraph(adj)  # Graph
    pr_ = page_rank(G, m2, n2)  # Page-Rank.
    LP = label_of_patch(edge_) 
    mask1, image=rankdown(im_copy, threshold)
    mask2, image=rankdown(-im_copy, threshold)
    Plot(pr_, edge_, im_copy, LP)
