import numpy as np
from itertools import permutations
from functools import reduce
import scipy.special
import multiprocessing
import scipy.optimize


def F(beta, S, alpha, gamma):
    return beta * S - alpha - gamma

def S(c, M):
    return (1 - c)**M

def cxtfun(hxtslice, r0, dx, n_x, Nh):
    # Number of x bins (antigenic space)
    B = n_x
    xvals = np.arange(0, B)
    fft_kernel = np.fft.fft(np.exp(-abs(xvals - int(B/2))* dx / r0))

    # Computing prob of infection through convolution
    fft_nh = np.fft.fft(hxtslice/(hxtslice.sum()*dx))
    cxt = np.fft.fftshift( np.real( np.fft.ifft( fft_nh*fft_kernel ) ) )*dx
    return cxt

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))/(np.sqrt(2*np.pi* np.power(sig,2.)))
    
def skewed_gaussian(x, sk, sigma, mu):
    """
    A skewed gaussian to initialize ntot
    """
    return np.exp(-(x-mu)**2/(sigma**2*2))*(1 +scipy.special.erf(sk*x/np.sqrt(2)))

def Jump_step(populations, jump_probas, i):
    
    return

































# END #