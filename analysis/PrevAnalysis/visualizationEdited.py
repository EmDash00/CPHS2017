import sys
import os
import glob
import importlib
#
import numpy as np
import pylab as plt
import matplotlib
import math as m
from matplotlib import rc
plt.close()
#
# Say, "the default sans-serif font is Helvetica"
rc('font',**{'sans-serif':'Helvetica','family':'sans-serif'})
rc('text',usetex=False)

help = """
usage:
  visualization subject protocol [formats]

add filename formats to save files in those formats:
  visualization subject protocol png,mp4
"""

def dbg(s):
  print s


# CHECK THAT ALL RELEVANT PARAMETERS ARE SPECIFIED IN COMMAND
args = sys.argv

if len(args) < 2:
  dbg('\nABORT -- no subject specified')
  dbg(help)
  sys.exit(0)

if len(args) < 3:
  dbg('\nABORT -- no protocol specified')
  dbg(help)
  sys.exit(0)

subject = args[1]
protocol = args[2]

fmts = []
if len(args) >= 4:
  fmts = args[3].split(',')

do_anim = True
do_anim = False
do_stills = True
#do_stills = False

# importing protocol (su17v3)
proto = importlib.import_module('protocols.'+protocol)

# find files for specific subject
fis = glob.glob(os.path.join('data','ref','*_'+protocol+'*'))
ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
ids = [id for id in ids if id[-3:] not in ['rej','oob','rst','max']]

from protocols.globals import ( size, RNG, SHIP_SHIFT, INP_SHIFT,
                                TIMES, _GOLD, _PURPLE,
                                XRNG, YRNG, GRID_SPACE )

#from util import *

INP_SHIFT = INP_SHIFT - SHIP_SHIFT
SHIP_SHIFT = SHIP_SHIFT - SHIP_SHIFT # 0.
XRNG = (-.5,2.5)

trials = dict()

for id in ids:
  #fis = sorted(glob.glob(os.path.join('data',subject,'*_'+protocol+'_'+id+'.npz')))
  fis = sorted(glob.glob(os.path.join('data','ref','*_'+protocol+'_'+id+'.npz')))
  if len(fis) > 1:
    dbg('WARNING -- repeated trials for id ='+id)
  assert len(fis) > 0, 'ERROR -- no data for id ='+id
  fi = fis[-1]
  print('LOAD '+fi)
  trial = dict(np.load(fi))
  trials[id] = trial

#sys.exit(0)
times_ = [trials[ids[0]]['time_']]
refs_  = [trials[ids[0]]['ref_']]
outs_  = [trials[ids[0]]['out_']]
inps_  = [trials[ids[0]]['inp_']]#/(scales[0]*3.)] # what's the scale for?
dists_ = [trials[ids[0]]['dis_']]

times = np.hstack(times_)[-2400:] # take out first 5 sec
refs = np.hstack(refs_)[-2400:]
outs = np.hstack(outs_)[-2400:]
inps = np.hstack(inps_)[-2400:]
dists = np.hstack(dists_)[-2400:]

samp0 = 60 # data being saved at 60 Hz = 1/60 sec/data point
timesplit = 1.0/samp0

# get frequencies
freq = np.linspace(0,599.5,refs.shape[0]/2)
tempfreq = np.ones((refs.shape[0],))
TEMPFREQ = np.fft.fft(tempfreq)[0]
OUTS = np.fft.fft(outs)[:outs.shape[0]/2]/TEMPFREQ
INPS = np.fft.fft(inps)[:inps.shape[0]/2]/TEMPFREQ
REFS = np.fft.fft(refs)[:refs.shape[0]/2]/TEMPFREQ
primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31])# max =37
IX = primes*2

Hru = np.divide(INPS[IX],REFS[IX])
Hru = np.asarray(Hru)

fis = glob.glob(os.path.join('data','dist','*_'+protocol+'*'))
ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
ids = [id for id in ids if id[-3:] not in ['rej','oob','rst','max']]

for id in ids:
  #fis = sorted(glob.glob(os.path.join('data',subject,'*_'+protocol+'_'+id+'.npz')))
  fis = sorted(glob.glob(os.path.join('data','dist','*_'+protocol+'_'+id+'.npz')))
  if len(fis) > 1:
    dbg('WARNING -- repeated trials for id ='+id)
  assert len(fis) > 0, 'ERROR -- no data for id ='+id
  fi = fis[-1]
  print('LOAD '+fi)
  trial = dict(np.load(fi))
  trials[id] = trial

#sys.exit(0)

times_ = [trials[ids[0]]['time_']]
refs_  = [trials[ids[0]]['ref_']]
outs_  = [trials[ids[0]]['out_']]
inps_  = [trials[ids[0]]['inp_']]#/(scales[0]*3.)] # what's the scale for?
dists_ = [trials[ids[0]]['dis_']]

times = np.hstack(times_)[-2400:]
refs = np.hstack(refs_)[-2400:]
outs = np.hstack(outs_)[-2400:]
inps = np.hstack(inps_)[-2400:]
dists = np.hstack(dists_)[-2400:]

# get frequencies
OUTS = np.fft.fft(outs)[:outs.shape[0]/2]/TEMPFREQ
INPS = np.fft.fft(inps)[:inps.shape[0]/2]/TEMPFREQ
REFS = np.fft.fft(refs)[:refs.shape[0]/2]/TEMPFREQ
DISTS = np.fft.fft(dists)[:dists.shape[0]/2]/TEMPFREQ

Hdu = np.divide(INPS[IX],DISTS[IX])
Hdu = np.asarray(Hdu)

M = 1/(1j*(2*m.pi*freq[IX]))*20
M = np.asarray(M)
B = np.divide(-Hdu,(np.multiply(M,(1.0+0j)+Hdu)))#np.divide((1+Hdu),np.multiply(M,Hdu))
B = np.asarray(B)

F = np.multiply(Hru,(1.0+0j)+np.multiply(B,M))-B#np.multiply(Hru,(np.multiply(np.ones(B.shape)+np.multiply(B,M)))))-B#np.divide(Hru,Hdu) - B
F= np.asarray(F)

# plot B and F

plt.figure(1)
plt.subplot(211)
plt.loglog(freq[IX],np.absolute(B),label='Back')
plt.loglog(freq[IX],np.absolute(F),label='Forward')
plt.ylabel('Gain')
plt.title('Estimated Feedback/Feedforward Policies')
plt.legend()

plt.subplot(212)
plt.semilogx(freq[IX],np.angle(B)*180/m.pi,label='Back')
plt.semilogx(freq[IX],np.angle(F)*180/m.pi,label='Forward')

#plot reference and output
# plt.figure(1)
# plt.subplot(241)
# plt.plot(times,refs,label='r')
# plt.plot(times,outs,label='y')
# plt.title('r and y time domain')
# plt.legend()
#
# # take fft
# plt.subplot(242)
# plt.plot(freq,np.absolute(REFS),label='R')
# plt.plot(freq,np.absolute(OUTS),label='Y')
# plt.title('R and Y freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.legend()
# plt.xlim(0,40)
# plt.subplot(246)
# plt.plot(freq[IX],np.angle(REFS[IX]),label='R')
# plt.plot(freq[IX],np.angle(OUTS[IX]),label='Y')
# plt.title('R and Y freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.legend()
# plt.xlim(0,40)
#
# plt.subplot(243)
# plt.plot(freq,np.absolute(INPS))
# plt.title('U freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(247)
# plt.plot(freq[IX],np.angle(INPS[IX]))
# plt.title('U freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)

# find H = Y/U = OUTS/INPS

# plt.subplot(244)
# plt.plot(freq[IX],np.absolute(Hru))
# plt.title('Hru magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(248)
# plt.plot(freq[IX],np.angle(Hru))
# #plt.plot(freq[0:2:70],np.angle(Hru[0:2:70]))
# plt.title('Hru phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)
#
#
#
# plt.figure(2)
#
# plt.subplot(241)
# plt.plot(times,refs,label='r')
# plt.plot(times,outs,label='y')
# plt.title('r and y time domain')
# plt.legend()
#
# # take fft
# plt.subplot(242)
# plt.semilogx(freq,20*np.log10(np.absolute(REFS)),label='R')
# plt.semilogx(freq,20*np.log10(np.absolute(OUTS)),label='Y')
# plt.title('R and Y freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.legend()
# plt.xlim(0,40)
# plt.subplot(246)
# plt.semilogx(freq[IX],np.angle(REFS[IX]),label='R')
# plt.semilogx(freq[IX],np.angle(OUTS[IX]),label='Y')
# plt.title('R and Y freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.legend()
# plt.xlim(0,40)
#
# plt.subplot(243)
# plt.semilogx(freq,20*np.log10(np.absolute(INPS)))
# plt.title('U freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(247)
# plt.semilogx(freq[IX],np.angle(INPS[IX]))
# plt.title('U freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)
#
# # find H = Y/U = OUTS/INPS
# plt.subplot(244)
# plt.semilogx(freq[IX],20*np.log10(np.absolute(Hru)))
# plt.title('Hru magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(248)
# plt.semilogx(freq[IX],np.angle(Hru))
# #plt.plot(freq[0:2:70],np.angle(Hru[0:2:70]))
# plt.title('Hru phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)

# find Hdu
# find files for specific subject

# plt.figure(3)
# plt.subplot(241)
# plt.plot(times,dists,label='d')
# plt.plot(times,outs,label='y')
# plt.title('d and y time domain')
# plt.legend()
#
# # take fft
# plt.subplot(242)
# plt.plot(freq,np.absolute(DISTS),label='D')
# plt.plot(freq,np.absolute(OUTS),label='Y')
# plt.title('D and Y freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.legend()
# plt.xlim(0,40)
# plt.subplot(246)
# plt.plot(freq[IX],np.angle(DISTS[IX]),label='D')
# plt.plot(freq[IX],np.angle(OUTS[IX]),label='Y')
# plt.title('D and Y freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.legend()
# plt.xlim(0,40)
#
# plt.subplot(243)
# plt.plot(freq,np.absolute(INPS))
# plt.title('U freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(247)
# plt.plot(freq[IX],np.angle(INPS[IX]))
# plt.title('U freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)

# find H = Y/U = OUTS/INPS

# plt.subplot(244)
# plt.plot(freq[IX],np.absolute(Hdu))
# plt.title('Hdu magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(248)
# plt.plot(freq[IX],np.angle(Hdu))
# #plt.plot(freq[0:2:70],np.angle(Hru[0:2:70]))
# plt.title('Hdu phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)
#
# plt.figure(4)
# plt.subplot(241)
# plt.plot(times,dists,label='d')
# plt.plot(times,outs,label='y')
# plt.title('d and y time domain')
# plt.legend()
#
# # take fft
# plt.subplot(242)
# plt.semilogx(freq,20*np.log10(np.absolute(DISTS)),label='D')
# plt.semilogx(freq,20*np.log10(np.absolute(OUTS)),label='Y')
# plt.title('R and Y freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.legend()
# plt.xlim(0,40)
# plt.subplot(246)
# plt.semilogx(freq[IX],np.angle(DISTS[IX]),label='D')
# plt.semilogx(freq[IX],np.angle(OUTS[IX]),label='Y')
# plt.title('D and Y freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.legend()
# plt.xlim(0,40)
#
# plt.subplot(243)
# plt.semilogx(freq,20*np.log10(np.absolute(INPS)))
# plt.title('U freq domain magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(247)
# plt.semilogx(freq[IX],np.angle(INPS[IX]))
# plt.title('U freq domain phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)
#
# # find H = Y/U = OUTS/INPS
# plt.subplot(244)
# plt.semilogx(freq[IX],20*np.log10(np.absolute(Hdu)))
# plt.title('Hdu magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
#
# plt.subplot(248)
# plt.semilogx(freq[IX],np.angle(Hdu))
# #plt.plot(freq[0:2:70],np.angle(Hru[0:2:70]))
# plt.title('Hdu phase')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
# plt.xlim(0,40)
#
# # FIND B - feedback section; B = (1+Hdu)/(MHdu)
# # create M
# plt.figure(5)

# plt.subplot(231)
# plt.plot(freq[IX],np.absolute(M))
# plt.title('M magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
# plt.subplot(234)
# plt.plot(freq[IX],np.angle(M))
# plt.title('M phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
#
# plt.subplot(232)
# plt.plot(freq[IX],np.absolute(B))
# plt.title('B magnitude')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.subplot(235)
# plt.plot(freq[IX],np.angle(B))
# plt.title('B phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')

# FIND F; F = Hru/Hdu - B

# plt.subplot(233)
# plt.plot(freq[IX],np.absolute(F))
# plt.title('F magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
# plt.subplot(236)
# plt.plot(freq[IX],np.angle(F))
# plt.title('F phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
#
#
# plt.figure(6)
# plt.subplot(231)
# plt.semilogx(freq[IX],20*np.log10(np.absolute(M)))
# plt.title('M magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude (dB)')
# plt.xlim(0,40)
# plt.subplot(234)
# plt.semilogx(freq[IX],np.angle(M))
# plt.title('M phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
#
# plt.subplot(232)
# plt.semilogx(freq[IX],20*np.log10(np.absolute(B)))
# plt.title('B magnitude')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.subplot(235)
# plt.semilogx(freq[IX],np.angle(B))
# plt.title('B phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
#
# # FIND F; F = Hru/Hdu - B
# plt.subplot(233)
# plt.semilogx(freq[IX],20*np.log10(np.absolute(F)))
# plt.title('F magnitude')
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Magnitude')
# plt.xlim(0,40)
# plt.subplot(236)
# plt.semilogx(freq[IX],np.angle(F))
# plt.title('F phase')
# plt.xlim(0,40)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase (radians)')
#
#
# # BOTH
# # find files for specific subject
# fis = glob.glob(os.path.join('data','both','*_'+protocol+'*'))
# ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
# ids = [id for id in ids if id[-3:] not in ['rej','oob','rst','max']]
#
# for id in ids:
#   #fis = sorted(glob.glob(os.path.join('data',subject,'*_'+protocol+'_'+id+'.npz')))
#   fis = sorted(glob.glob(os.path.join('data','both','*_'+protocol+'_'+id+'.npz')))
#   if len(fis) > 1:
#     dbg('WARNING -- repeated trials for id ='+id)
#   assert len(fis) > 0, 'ERROR -- no data for id ='+id
#   fi = fis[-1]
#   print('LOAD '+fi)
#   trial = dict(np.load(fi))
#   trials[id] = trial
#
# #sys.exit(0)
#
# times_ = [trials[ids[0]]['time_']]
# refs_  = [trials[ids[0]]['ref_']]
# outs_  = [trials[ids[0]]['out_']]
# inps_  = [trials[ids[0]]['inp_']]#/(scales[0]*3.)] # what's the scale for?
# dists_ = [trials[ids[0]]['dis_']]
#
# times = np.hstack(times_)[-2400:]
# refs = np.hstack(refs_)[-2400:]
# outs = np.hstack(outs_)[-2400:]
# inps = np.hstack(inps_)[-2400:]
# dists = np.hstack(dists_)[-2400:]
#
# samp0 = 60 # data being saved at 60 Hz = 1/60 sec/data point
# timesplit = 1.0/samp0
#
# freq = np.linspace(0,599.5,refs.shape[0]/2)
# tempfreq = np.ones((refs.shape[0],))
# TEMPFREQ = np.fft.fft(tempfreq)[0]
# primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31])# max =37
# primesE = np.asarray([2,5,11,17,23,31])#primes[::2]
# primesO = np.asarray([3,7,13,19,29])#primes[1::2]
# IXE = primesE*2
# IXO = primesO*2
# IX = primes*2
# # REFS - E, DIS - O
# OUTS1 = np.fft.fft(outs)[:outs.shape[0]/2]/TEMPFREQ
# INPS1 = np.fft.fft(inps)[:inps.shape[0]/2]/TEMPFREQ
# REFS1 = np.fft.fft(refs)[:refs.shape[0]/2]/TEMPFREQ
# DISTS1 = np.fft.fft(dists)[:dists.shape[0]/2]/TEMPFREQ
#
# times_ = [trials[ids[1]]['time_']]
# refs_  = [trials[ids[1]]['ref_']]
# outs_  = [trials[ids[1]]['out_']]
# inps_  = [trials[ids[1]]['inp_']]#/(scales[0]*3.)] # what's the scale for?
# dists_ = [trials[ids[1]]['dis_']]
#
# times = np.hstack(times_)[-2400:]
# refs = np.hstack(refs_)[-2400:]
# outs = np.hstack(outs_)[-2400:]
# inps = np.hstack(inps_)[-2400:]
# dists = np.hstack(dists_)[-2400:]
# # REFS - O, DIS - E
# OUTS2 = np.fft.fft(outs)[:outs.shape[0]/2]/TEMPFREQ
# INPS2 = np.fft.fft(inps)[:inps.shape[0]/2]/TEMPFREQ
# REFS2 = np.fft.fft(refs)[:refs.shape[0]/2]/TEMPFREQ
# DISTS2 = np.fft.fft(dists)[:dists.shape[0]/2]/TEMPFREQ
#
# #print np.concatenate([REFS1[IXE],REFS2[IXO]])
# REFS = []
# print REFS2[IX]
# print DISTS2[IX]
# for kk in range(len(IXO)):
#     #temp = [REFS1[IXE][kk],REFS2[IXO][kk]]
#     REFS.append(REFS1[IXE][kk])
#     REFS.append(REFS2[IXO][kk])
#     #REFS.append(REFS2[IXO][kk])
# REFS.append(REFS1[IXE][-1])
# REFS = np.hstack(REFS)
# #print REFS
# #REFS = np.concatenate([REFS1[IXE],REFS2[IXO]])
# OUTS = np.concatenate([OUTS1[IXE],OUTS2[IXO]])
# INPS = np.concatenate([INPS1[IXE],INPS2[IXO]])
# DISTS = np.concatenate([DISTS1[IXE],DISTS2[IXO]])
# freqAll = np.concatenate([IXE,IXO])
# # find refs
# # M = 1/(1j*(2*m.pi*freq[IX]))
# # Hru = np.divide(INPS[IX],REFS[IX])
# # Bboth = np.divide((1+Hdu),np.multiply(M,Hdu))
# # Fboth = np.divide(Hru,Hdu) - Bboth
#
# #plot reference and output
# plt.figure(7)
#
# # plt.subplot(241)
# # plt.plot(times,refs,label='r')
# # plt.plot(times,outs,label='y')
# # plt.title('r and y time domain')
# # plt.legend()
#
# # take fft
# plt.subplot(242)
# plt.plot(freq[IX],np.absolute(REFS),label='R')
# # plt.plot(freq,np.absolute(OUTS),label='Y')
# # plt.title('R and Y freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.legend()
# # plt.xlim(0,40)
# # plt.subplot(246)
# # plt.plot(freq[IX],np.angle(REFS[IX]),label='R')
# # plt.plot(freq[IX],np.angle(OUTS[IX]),label='Y')
# # plt.title('R and Y freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.legend()
# # plt.xlim(0,40)
# #
# # plt.subplot(243)
# # plt.plot(freq,np.absolute(INPS))
# # plt.title('U freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(247)
# # plt.plot(freq[IX],np.angle(INPS[IX]))
# # plt.title('U freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # # find H = Y/U = OUTS/INPS
# # plt.subplot(244)
# # plt.plot(freq[IX],np.absolute(Hru))
# # plt.title('Hru magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(248)
# # plt.plot(freq[IX],np.angle(Hru))
# # plt.title('Hru phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# #
# #
# # plt.figure(8)
# #
# # plt.subplot(241)
# # plt.plot(times,refs,label='r')
# # plt.plot(times,outs,label='y')
# # plt.title('r and y time domain')
# # plt.legend()
# #
# # # take fft
# # plt.subplot(242)
# # plt.semilogx(freq,20*np.log10(np.absolute(REFS)),label='R')
# # plt.semilogx(freq,20*np.log10(np.absolute(OUTS)),label='Y')
# # plt.title('R and Y freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.legend()
# # plt.xlim(0,40)
# # plt.subplot(246)
# # plt.semilogx(freq[IX],np.angle(REFS[IX]),label='R')
# # plt.semilogx(freq[IX],np.angle(OUTS[IX]),label='Y')
# # plt.title('R and Y freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.legend()
# # plt.xlim(0,40)
# #
# # plt.subplot(243)
# # plt.semilogx(freq,20*np.log10(np.absolute(INPS)))
# # plt.title('U freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(247)
# # plt.semilogx(freq[IX],np.angle(INPS[IX]))
# # plt.title('U freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # # find H = Y/U = OUTS/INPS
# # Hru = np.divide(INPS[IX],REFS[IX])
# # plt.subplot(244)
# # plt.semilogx(freq[IX],20*np.log10(np.absolute(Hru)))
# # plt.title('Hru magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(248)
# # plt.semilogx(freq[IX],np.angle(Hru))
# # #plt.plot(freq[0:2:70],np.angle(Hru[0:2:70]))
# # plt.title('Hru phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # plt.figure(9)
# # plt.subplot(241)
# # plt.plot(times,dists,label='d')
# # plt.plot(times,outs,label='y')
# # plt.title('d and y time domain')
# # plt.legend()
# #
# # # take fft
# # plt.subplot(242)
# # plt.plot(freq,np.absolute(DISTS),label='D')
# # plt.plot(freq,np.absolute(OUTS),label='Y')
# # plt.title('D and Y freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.legend()
# # plt.xlim(0,40)
# # plt.subplot(246)
# # plt.plot(freq[IX],np.angle(DISTS[IX]),label='D')
# # plt.plot(freq[IX],np.angle(OUTS[IX]),label='Y')
# # plt.title('D and Y freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.legend()
# # plt.xlim(0,40)
# #
# # plt.subplot(243)
# # plt.plot(freq,np.absolute(INPS))
# # plt.title('U freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(247)
# # plt.plot(freq[IX],np.angle(INPS[IX]))
# # plt.title('U freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # # find H = Y/U = OUTS/INPS
# # plt.subplot(244)
# # plt.plot(freq[IX],np.absolute(Hdu))
# # plt.title('Hdu magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(248)
# # plt.plot(freq[IX],np.angle(Hdu))
# # plt.title('Hdu phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # plt.figure(10)
# # plt.subplot(241)
# # plt.plot(times,dists,label='d')
# # plt.plot(times,outs,label='y')
# # plt.title('d and y time domain')
# # plt.legend()
# #
# # # take fft
# # plt.subplot(242)
# # plt.semilogx(freq,20*np.log10(np.absolute(DISTS)),label='D')
# # plt.semilogx(freq,20*np.log10(np.absolute(OUTS)),label='Y')
# # plt.title('R and Y freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.legend()
# # plt.xlim(0,40)
# # plt.subplot(246)
# # plt.semilogx(freq[IX],np.angle(DISTS[IX]),label='D')
# # plt.semilogx(freq[IX],np.angle(OUTS[IX]),label='Y')
# # plt.title('D and Y freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.legend()
# # plt.xlim(0,40)
# #
# # plt.subplot(243)
# # plt.semilogx(freq,20*np.log10(np.absolute(INPS)))
# # plt.title('U freq domain magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(247)
# # plt.semilogx(freq[IX],np.angle(INPS[IX]))
# # plt.title('U freq domain phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# # # find H = Y/U = OUTS/INPS
# # plt.subplot(244)
# # plt.semilogx(freq[IX],20*np.log10(np.absolute(Hdu)))
# # plt.title('Hdu magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# #
# # plt.subplot(248)
# # plt.semilogx(freq[IX],np.angle(Hdu))
# # plt.title('Hdu phase')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# # plt.xlim(0,40)
# #
# #
# # plt.figure(11)
# # plt.subplot(231)
# # plt.plot(freq[IX],np.absolute(M))
# # plt.title('M magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# # plt.subplot(234)
# # plt.plot(freq[IX],np.angle(M))
# # plt.title('M phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# #
# # plt.subplot(232)
# # plt.plot(freq[IX],np.absolute(Bboth))
# # plt.title('B magnitude both')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.subplot(235)
# # plt.plot(freq[IX],np.angle(Bboth))
# # plt.title('B both phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# #
# # # FIND F; F = Hru/Hdu - B
# # plt.subplot(233)
# # plt.plot(freq[IX],np.absolute(Fboth))
# # plt.title('F magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# # plt.subplot(236)
# # plt.plot(freq[IX],np.angle(Fboth))
# # plt.title('F phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
#
#
#
#
#
#
#
#
#
#
# #
# # plt.figure(6)
# # plt.subplot(231)
# # plt.semilogx(freq[IX],20*np.log10(np.absolute(M)))
# # plt.title('M magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude (dB)')
# # plt.xlim(0,40)
# # plt.subplot(234)
# # plt.semilogx(freq[IX],np.angle(M))
# # plt.title('M phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# #
# # plt.subplot(232)
# # plt.semilogx(freq[IX],20*np.log10(np.absolute(Bboth)))
# # plt.title('B magnitude')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.subplot(235)
# # plt.semilogx(freq[IX],np.angle(Bboth))
# # plt.title('B phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
# #
# # # FIND F; F = Hru/Hdu - B
# # F = np.divide(Hru,Hdu) - B
# # plt.subplot(233)
# # plt.semilogx(freq[IX],20*np.log10(np.absolute(Fboth)))
# # plt.title('F magnitude')
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0,40)
# # plt.subplot(236)
# # plt.semilogx(freq[IX],np.angle(Fboth))
# # plt.title('F phase')
# # plt.xlim(0,40)
# # plt.xlabel('Freq (Hz)')
# # plt.ylabel('Phase (radians)')
#
# # #
# # # plt.figure(5)
# # # plt.subplot(221)
# # # plt.semilogx(freq,20*np.log10(np.absolute(Bboth)))#,marker='o')
# # # plt.title('B magnitude')
# # # #plt.xlim(0,20)
# # # plt.subplot(223)
# # # plt.semilogx(freq,np.angle(Bboth))
# # # plt.title('B phase')
# # #
# # #
# # # # FIND F; F = Hru/Hdu - B
# # # plt.subplot(222)
# # # plt.semilogx(freq,20*np.log10(np.absolute(Fboth)))
# # # plt.title('F magnitude')
# # # plt.subplot(224)
# # # plt.semilogx(freq,np.angle(Fboth))
# # # plt.title('F phase')
# # #
# # #
# # # # similarity between F and B
# # # Fdiffmag = (np.absolute(Fboth) - np.absolute(F))
# # # Fdiffphase = (np.angle(Fboth) - np.angle(F))
# # #
# # # Bdiffmag = (np.absolute(Bboth) - np.absolute(B))
# # # Bdiffphase = (np.angle(Bboth) - np.angle(B))
# # #
# # # plt.figure(6)
# # # plt.subplot(221)
# # # plt.semilogx(freq,20*np.log10(Bdiffmag))#,marker='o')
# # # plt.title('B diff magnitude')
# # # #plt.xlim(0,20)
# # # plt.subplot(223)
# # # plt.semilogx(freq,Bdiffphase)
# # # plt.title('B diff phase')
# # #
# # #
# # # # FIND F; F = Hru/Hdu - B
# # # plt.subplot(222)
# # # plt.semilogx(freq,20*np.log10(Fdiffmag))
# # # plt.title('F diff magnitude')
# # # plt.subplot(224)
# # # plt.semilogx(freq,Fdiffphase)
# # # plt.title('F diff phase')

plt.show()
