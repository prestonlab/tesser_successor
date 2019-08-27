# creating and testing model RDMs
import numpy as np
import scipy.spatial.distance as spsd
import matplotlib.pyplot as plt
import util
import sr


sr_data, dl, induc_data, idl = util.read_file()

# Minor function to remove zeros in matrices with missing objects
def fix(matrix):
    new_matrix = matrix
    new_matrix[new_matrix == 0] = 0.0000001
    return new_matrix

# Transposing a matrix pattern to check for symmetry
def trans(matrix):
    m = np.transpose(matrix)
    return m

#
def make_rdm (matrix):
    return spsd.squareform (spsd.pdist (matrix, 'correlation'))

#
def track_alpha_with_learning_rmd_change(subj_num,gamma,alpha_min ,alpha_max ,step):
    for x in np.arange(alpha_min, alpha_max, step):
        gamma = np.round(gamma,decimals=2)
        alpha = np.round(x,decimals=2)
        r = 7.5
        c = 5
        M = np.zeros([21,21])  
        a = util.initial_subj_run(subj_num)
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(r, c))
        plt.suptitle("Trancking RDM with Learning and " r"$\alpha$ "" Change " + "\n"
                     + "Subject Number: " + str(subj_num) +" Subject Name: " + dl[a][:17] + "  "r'$\alpha$ : ' + str(alpha) + "\n", 
                 fontweight="bold", size = (10), ha = 'center',y=1.10)
        plt.subplots_adjust(hspace=1)
        title = " Part_" + str(1)
        for run in [0,1,2,3,4]:
            if run < 3:
                row = 0
                ax = axes[row,run]
                M_new = np.copy(M)
                M_new = fix(sr.make_M_Matrix((a + run), gamma, alpha, np.copy(M)))
                rdm = make_rdm(M_new)
                ax.set_title(title + " Run_" + dl[a + run][-5] + "\n")
                im = ax.imshow(rdm, cmap='viridis', vmin=0, vmax=1)
                M = M_new

            if run > 2:
                row = 1
                ax = axes[row,run-3]
                M_new = np.copy(M)
                M_new = fix(sr.make_M_Matrix((a + run), gamma, alpha, np.copy(M)))
                rdm = make_rdm(M_new)
                ax.set_title(title + " Run_" + dl[a + run][-5] + "\n")
                im = ax.imshow(rdm, cmap='viridis', vmin=0, vmax=1)
                M = M_new
                
            if run == 4:
                fig.delaxes(axes[1][2])

        cbar = fig.colorbar(im, ax=axes.ravel().tolist())
        cbar.set_label('Color Bar')
       
    return plt.show()