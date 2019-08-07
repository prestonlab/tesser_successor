import numpy as np
import matplotlib.pyplot as plt
import sr


def plot_explore_runs(OPTION, SUBJECT, GAMMA, ALPHA):
    figsize=(70,30)
    size = 50    
    fig, axs = plt.subplots(2, 6, figsize=figsize ,sharex='col', sharey='row')
    plt.suptitle('Learning: ' + OPTION +'  SUBJECT: ' + str(SUBJECT) + " with "r'$\gamma$ : ' + str(GAMMA) + " and "r'$\alpha$ : ' + str(ALPHA),size = 80)
    M,pr = sr.explore_runs(OPTION, SUBJECT, GAMMA, ALPHA)
    for i in range(len(M)):
        part = int(pr[i][0])
        run = int(pr[i][1])
        axs[part - 1, run-1].matshow(M[i], cmap='viridis', vmin=0, vmax=1)
        axs[part - 1, run -1].set_title('Part_%s Run_%s'% (pr[i][0],pr[i][1]), size=size)


    fig.delaxes(axs[0][5])
    return plt.show()