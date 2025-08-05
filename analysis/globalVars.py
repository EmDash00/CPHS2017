import numpy as np
import math as m

# Define things
primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19])# max =37
IX = primes*2
primesE = np.asarray([2,5,11,17])#primes[::2]
primesO = np.asarray([3,7,13,19])#primes[1::2]
IXE = primesE*2
IXO = primesO*2
freq = np.linspace(0,599.5,1200)#refs.shape[0]/2)
tempfreq = np.ones((1200,))
TEMPFREQ = np.fft.fft(tempfreq)[0]
s = 1.0j*2*m.pi*freq[IX]/20
M = 1.0/(np.multiply(s,s)+s) # plus or minus s??? # TODO
M = np.asarray(M)
