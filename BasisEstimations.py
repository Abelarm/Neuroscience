__author__ = 'Luigi'


import hrf_estimation
import numpy as np
from scipy.io import savemat

def main():

    t = np.array([0,0,1,2,3,4,5,6], int)


    #print t.shape

    #basis1=np.arange(8)
    #print basis1.shape

    basis1=np.array(hrf_estimation.hrf.spmt(t),float)
    basis2=np.array(hrf_estimation.hrf.dspmt(t),float)
    basis3=np.array(hrf_estimation.hrf.ddspmt(t),float)

    toreturn=np.array([basis1,basis2,basis3],float)


    savemat('basis.mat',mdict={'basis':toreturn})


if __name__ == '__main__':
    main()