import matplotlib.pyplot as plt
import numpy as np

from seqer import *

#WWWWWWWWDDDDDDDDDRRRRRRRRRKK
#plt.style.use('_mpl-gallery')

def plot_charges(sequence):
    # make the data
    x = np.arange(start=1, stop=len(str(sequence))+1, step=1)

    y = []
    for i in range(len(str(sequence))):
        y.append(sequence.indiv_charge(i))

    
    plt.plot(x, y)
    plt.xlabel('Residue Placement')
    plt.ylabel('Charge')

    plt.show()
