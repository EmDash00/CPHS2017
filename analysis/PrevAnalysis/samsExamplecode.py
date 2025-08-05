# from _.py

import sys
import os
import glob
import importlib
from copy import deepcopy
#
import numpy as np
from scipy import linalg as la
import pylab as plt
import matplotlib
from matplotlib import rc
#
# Say, "the default sans-serif font is Helvetica"
rc('font',**{'sans-serif':'Helvetica','family':'sans-serif'})
rc('text',usetex=False)
#
# Say, "only print 2 digits after decimal"
np.set_printoptions(precision=2)

if 1:

  from references import sum_of_sines_ramp as sos
  from references import zero as zero
  refs = dict(sos=sos,zer=zero)

  # vector fields
  #vfs = ['so']
  vfs = ['fo']
  #vfs = ['fo','so']
  #
  states = dict(fo=[0.],so=[0.,0.])
  # time to complete task
  dt = 1.
  #seed = 49
  #np.random.seed(seed)
  #
  period = 20 # sec
  f_base = 1./period # Hz
  # TODO go up to 3--5Hz -- check that Sam can track, and check when becomes sub-pixel
  primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199])
  f_max = dict(fo=2,so=1)
  p_max = dict([(vf,np.nonzero(primes*f_base <= f_max[vf])[0][-1]) for vf in vfs])
  f_primes = dict([(vf,np.asarray(primes[:p_max[vf]])) for vf in vfs])
  #
  frequencies = dict([(vf,f_primes[vf]*f_base) for vf in vfs])
  frequencies_r = frequencies.copy()
  frequencies_d = frequencies.copy()
  #
  #print 'p_max = ',p_max,', primes[p_max] = ',primes[p_max]
  amplitudes = dict([(vf,1./f_primes[vf]/f_primes[vf].size) for vf in vfs])
  amplitudes = dict([(vf,(1./f_primes[vf])*(0.5/f_primes[vf]).sum()) for vf in vfs])
  amplitudes_r = amplitudes.copy()
  amplitudes_d = amplitudes.copy()
  #
  num_refs = 5
  phase_shifts_r = dict([(vf,np.random.rand(num_refs,p_max[vf])) for vf in vfs])
  # only first disturbance signal differs from reference signals
  phase_shifts_d = deepcopy(phase_shifts_r)
  for vf in vfs:
    phase_shifts_d[vf][0] = np.random.rand(p_max[vf])
  #
  ramp = 0.25*period
  duration = 2*period + ramp
  #
  trial = dict(ramp=ramp,duration=duration)

  plt.figure(1); plt.clf()
  t = np.arange(ramp,ramp+period,period/(2.**10))

  vf = 'fo'
  sines_r = np.arange(p_max[vf])
  sines_d = np.arange(p_max[vf])

  ref_ = 'sos'
  dis_ = 'sos'

  for shift_id in range(num_refs):
    ref = lambda t,trial,x=None : refs[ref_[:3]](t,trial,x,
                                             frequencies_r[vf][sines_r],
                                             amplitudes_r[vf][sines_r],
                                             phase_shifts_r[vf][shift_id][sines_r])
    dis = lambda t,trial,x=None : refs[dis_[:3]](t,trial,x,
                                             frequencies_d[vf][sines_d],
                                             amplitudes_d[vf][sines_d],
                                             phase_shifts_d[vf][shift_id][sines_d])
    #
    trial.update(dict(vf=vf,ref=ref,dis=dis))
    #
    r = ref(t,trial)
    if np.any(np.abs(r) > 1.0):
      print "LARGE"
    plt.plot(t,r)

  plt.ylim(-1,+1)

  sys.exit(0)


if 0:
  seed = 49
  #np.random.seed(seed)
  #
  period = 20 # sec
  f_base = 1./period # Hz
  # TODO go up to 3--5Hz -- check that Sam can track, and check when becomes sub-pixel
  primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199])
  f_max = 2 # max Hz of signal
  p_max = np.nonzero(primes*f_base <= f_max)[0][-1]
  f_primes = np.asarray(primes[:p_max]) # prime multiples of f_base frequency
  print 'p_max = ',p_max,', primes[p_max] = ',primes[p_max]
  a = 1./f_primes;
  a = a / f_primes.size
  #a /= a.sum() # confine signal to [-1,+1]
  # TODO choose 3 random phases, apply negative of each
  ramp = 0.
  duration = 1*period

  plt.figure(1); plt.clf()
  t = np.arange(ramp,ramp+period,period/(2.**10))
  from protocols.references import sum_of_sines as sos

  for _ in range(5):
    p = np.random.rand(len(f_primes)) # phase ranges from 0 to 1
    trial = dict(frequencies=f_primes*f_base,amplitudes=a,phase_shifts=p,
                 ramp=ramp,duration=duration)
    r = sos(t,trial)
    if np.any(np.abs(r) > 1.0):
      print "LARGE"
    plt.plot(t,r)

  plt.ylim(-1,+1)

  sys.exit(0)

if 1:

  #import control as ctrl

  def sim(num,den,inp,num_=[]):
    """
    simulate SISO transfer function

    inputs:
      num - 1 x o_num - numerator z-transform coeff's descending to order zero
      den - 1 x o_den - denominator z-transform coeff's in descending order
      inp - 1 x n - input
      (optional)
      num_ - 1 x o_num_ - acausal numerator z-transform coeff's descending from order zero

    outputs:
      out - 1 x n - output:  out = (num[+ num_] / den) inp
    """
    o_num = len(num)
    o_den = len(den)
    o = max(o_num,o_den)
    o_num_ = len(num_)
    n = len(inp)
    #
    num = np.hstack((num,num_))
    den = np.asarray(den)
    #
    inp = np.hstack([np.zeros(o-1),inp,np.zeros(o_num_)])
    out = np.zeros_like(inp)
    #
    #1/0
    for k in range(o,n+o):
      out[k-1] = ( (num * inp[k-o_num:k+o_num_]).sum()
                 - (den[:-1] * out[k-o_den:k-1]).sum() ) / den[-1]

    #return out[o-1:]
    return inp,out

  def fit(inp,out,o_num,o_den,o_num_=0):
    #
    o = max(o_num,o_den)
    #
    h_inp = la.hankel(inp)[:-o+1,o-o_num:o+o_num_]
    h_out = la.hankel(out)[:-o+1,o-o_den:o]
    #
    #A = np.hstack((h_inp,h_out[:,:-1]))
    #b = h_out[:,-1:]
    #x = np.dot( la.inv(np.dot(A.T,A)), np.dot(A.T,b) ).flatten()
    #
    A = np.hstack((h_inp,-h_out[:,:-1]))
    b = h_out[:,-1:]
    x = la.lstsq(A,b)[0].flatten()
    #
    if o_num_ == 0:
      num,den = x[:o_num],np.hstack((x[o_num:],1.))
      return num,[],den
    #
    else:
      num,num_,den = x[:o_num],x[o_num:o_num+o_num_],np.hstack((x[o_num+o_num_:],1.))
      return num,num_,den


  o_num = 2
  o_den = 3
  o_num_ = 0

  num = np.random.randn(o_num)
  if o_num > 1:
    num = np.poly(np.roots(num[::-1])/(1e-0+np.abs(np.roots(num[::-1])).max()))[::-1]
  den = np.random.randn(o_den)
  den = np.poly(np.roots(den[::-1])/(1e-0+np.abs(np.roots(den[::-1])).max()))[::-1]
  num /= den[-1]
  den /= den[-1]
  num_ = np.random.randn(o_num_)

  assert np.all(np.abs(np.roots(num[::-1])) < 1)
  assert np.all(np.abs(np.roots(den[::-1])) < 1)

  F_ = 1./10
  T = 2./F_
  t = np.arange(0,T,1./20)
  primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199])
  F = F_*primes[primes < 30]
  A = 1./(2*np.pi*F)
  P = 2*np.pi*np.random.rand(F.size)

  from protocols.references import sum_of_sines as sos

  inp_ = sos(t,None,None,F,A,P)

  inp,out = sim(num,den,inp_,num_=num_)

  noi = 0e-1*np.random.randn(out.size)

  o = max(o_num,o_den)
  #
  #h_inp = la.hankel(inp)[:-o+1,o-o_num:o+o_num_]
  #h_out = la.hankel(out)[:-o+1,o-o_den:o]
  #
  #A = np.hstack((h_inp,h_out[:,:-1]))
  #b = h_out[:,-1:]
  ##x = np.dot( la.inv(np.dot(a.T,a)), np.dot(a.T,b) ).flatten()
  #x = la.lstsq(A,b)[0].flatten()
  #_num,_num_,_den = x[:o_num],x[o_num:o_num+o_num_],np.hstack((x[o_num+o_num_:],1.))

  _num,_num_,_den = fit(inp,out+noi,o_num,o_den,o_num_)
  _num,_num_,_den = fit(inp,out+noi,5,5,o_num_)
  _inp,_out = sim(_num,_den,inp_,num_=_num_)

  plt.figure(1); plt.clf()
  ax = plt.subplot(1,1,1); ax.grid('on')
  ax.plot(t,out[-len(t):],'g-',lw=1)
  ax.plot(t,(out+noi)[-len(t):],'g.')
  ax.plot(t,_out[-len(t):],'g--',lw=3)
  #ax.set_ylabel('output')
  #ax = plt.subplot(2,1,2); ax.grid('on')
  ax.plot(t,inp[-len(t):],'-',lw=1)
  #ax.plot(t,_inp[-len(t):],'--',lw=4)
  #ax.set_ylabel('input')
  ax.set_xlabel('time')
  plt.tight_layout()

  print 'CLEAN'
  #print ('allclose? ',
  #       np.allclose(num / den[-1], _num),
  #       np.allclose(num_ / den[-1], _num_),
  #       np.allclose(den / den[-1], _den))
  print
  print 'num'
  #print num / den[-1]
  #print _num
  print np.roots((num / den[-1])[::-1])
  print np.roots((_num)[::-1])
  print
  print 'num_'
  #print num_ / den[-1]
  #print _num_
  print np.roots((num_ / den[-1])[::-1])
  print np.roots((_num_)[::-1])
  print
  print 'den'
  #print den / den[-1]
  #print _den
  print np.roots((den / den[-1])[::-1])
  print np.roots((_den)[::-1])
  print

  sys.exit(0)

  for _ in range(2):
    inp_ = inp + 1e-4*np.random.randn(inp.size)
    out_ = out + 1e-4*np.random.randn(out.size)

    h_inp_ = la.hankel(inp_)[:-o+1,o-o_num:o]
    h_out_ = la.hankel(out_)[:-o+1,o-o_den:o]

    a = np.hstack((h_inp_,h_out_[:,:-1]))
    b = h_out_[:,-1:]
    x = np.dot( la.inv(np.dot(a.T,a)), np.dot(a.T,b) ).flatten()
    num__,den__ = x[:o_num],np.hstack((x[o_num:],1.))

    print 'NOISY'
    print 'num'
    print num / den[-1]
    print num__
    print
    print 'den'
    print den / den[-1]
    print den__
    print

  sys.exit(0)

if 1:
  d = dict(np.load('data/sam-su17/20170720-102207_su17_fo_r-zero_d-sos_s49_i0.npz'))

  def fit(I,x,y):
    """
    solves x a = y indexing into x via i:
    """
    I = [np.asarray(i) for i in I]
    x = np.copy(x); y = np.copy(y).flatten()
    m = min(0,np.hstack(I).min()); M = max(0,np.hstack(I).max())
    if len(x.shape) == 1:
      x.shape = (1,x.size)
    elif len(x.shape) == 2 and x.shape[0] == y.size and not x.shape[1] == y.size:
      x = x.T
    #assert x.shape[1] >= y.size
    #
    N = y.size
    X = []
    Y = []
    #
    for n in range(N):
      if n + m < 0 or n + M >= N:
        continue
      Xn = []
      for r,i in enumerate(I):
        Xn.append(x[r,n+i])
      X.append(np.hstack(Xn))
      Y.append([y[n]])
    #
    X = np.asarray(X)
    Y = np.asarray(Y)
    #
    a = np.dot(la.inv(np.dot(X.T,X)),np.dot(X.T,Y))
    #1/0
    return a.flatten()

  N = 50
  #x = np.sin(np.random.rand(N)*2*np.pi)
  x = np.sin(np.linspace(0.,3*np.pi/2,N))

  n = 3
  mode = 'same'
  mode = 'valid'
  for a in [[-1,2]#[-1,1],
            #np.random.randn(n),
            #[1,-2,1],
            ]:
    #I = [np.arange(-len(a),0)+1]
    I = [np.arange(-10,0)+1]
    y = np.convolve(x,a[::-1],mode)
    #a_ = fit(I,x,np.hstack([0.,y]))
    #y_ = np.convolve(x,a_[::-1],mode)
    X = np.asarray([x[r+I[0]] for r in range(1,N)])
    Y = y
    a__ = np.dot(la.inv(np.dot(X.T,X)),np.dot(X.T,Y))
    y__ = np.convolve(x,a__[::-1],mode)
    print 'convolve:'
    print a
    #print a_
    print a__
    #print 'a_ allclose? ', np.allclose(a, a_,rtol=1e-1,atol=1e-2)
    #print 'a__ allclose? ',np.allclose(a, a__,rtol=1e-1,atol=1e-2)
    #print 'y_ allclose? ', np.allclose(y, y_,rtol=1e-1,atol=1e-2)
    #print 'y__ allclose? ',np.allclose(y, y__,rtol=1e-1,atol=1e-2)
    #print y
    #print y_

  sys.exit(0)
