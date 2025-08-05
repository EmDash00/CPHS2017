#!/usr/bin/python
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
from scipy.signal import coherence

# define some variables that we need
global protocol
freq = np.linspace(0,599.5,1200)
tempfreq = np.ones((2400,))
TEMPFREQ = np.fft.fft(tempfreq)[0]
protocol = 'su17v3'

def findFFT(listofIDs,IX,subject):
    trials = dict()
    for id in listofIDs:
      fis = sorted(glob.glob(os.path.join('data',subject,'*_'+protocol+'_'+id+'.npz')))
      if len(fis) > 1:
        dbg('WARNING -- repeated trials for id ='+id)
      assert len(fis) > 0, 'ERROR -- no data for id ='+id
      fi = fis[-1]
      print('LOAD '+fi)
      trial = dict(np.load(fi))
      trials[id] = trial

    timedomainvalues = {}
    timedomainvalues[listofIDs[0]] = {}
    times_ = [trials[listofIDs[0]]['time_']]
    refs_  = [trials[listofIDs[0]]['ref_']]
    outs_  = [trials[listofIDs[0]]['out_']]
    inps_  = [trials[listofIDs[0]]['inp_']]#/(scales[0]*3.)] # what's the scale for?
    dists_ = [trials[listofIDs[0]]['dis_']]


    timedomainvalues[listofIDs[0]]['times'] = np.hstack(times_)[-2400:] # take out first 5 sec
    timedomainvalues[listofIDs[0]]['refs'] = np.hstack(refs_)[-2400:]
    timedomainvalues[listofIDs[0]]['outs'] = np.hstack(outs_)[-2400:]
    timedomainvalues[listofIDs[0]]['inps'] = np.hstack(inps_)[-2400:]
    timedomainvalues[listofIDs[0]]['dists'] = np.hstack(dists_)[-2400:]

    for id in listofIDs[1:]:
      times_ = [trials[id]['time_']]
      refs_  = [trials[id]['ref_']]
      outs_  = [trials[id]['out_']]
      inps_  = [trials[id]['inp_']]
      dists_ = [trials[id]['dis_']]
      timedomainvalues[id] = {}
      timedomainvalues[id]['times'] = (np.hstack(times_)[-2400:]) # take out first 5 sec
      timedomainvalues[id]['refs']=(np.hstack(refs_)[-2400:])
      timedomainvalues[id]['outs']=(np.hstack(outs_)[-2400:])
      timedomainvalues[id]['inps']=(np.hstack(inps_)[-2400:])
      timedomainvalues[id]['dists']=(np.hstack(dists_)[-2400:])


    FREQDOMAINVALUES = {}
    for id in listofIDs:
        FREQDOMAINVALUES[id] = {}
        FREQDOMAINVALUES[id]['OUTS'] = np.fft.fft(timedomainvalues[id]['outs'])[:timedomainvalues[id]['outs'].shape[0]/2]/TEMPFREQ
        FREQDOMAINVALUES[id]['INPS'] = np.fft.fft(timedomainvalues[id]['inps'])[:timedomainvalues[id]['inps'].shape[0]/2]/TEMPFREQ
        FREQDOMAINVALUES[id]['REFS'] = np.fft.fft(timedomainvalues[id]['refs'])[:timedomainvalues[id]['refs'].shape[0]/2]/TEMPFREQ
        FREQDOMAINVALUES[id]['DISTS'] = np.fft.fft(timedomainvalues[id]['dists'])[:timedomainvalues[id]['dists'].shape[0]/2]/TEMPFREQ
        FREQDOMAINVALUES[id]['Hru'] = np.divide(FREQDOMAINVALUES[id]['INPS'][IX],FREQDOMAINVALUES[id]['REFS'][IX])
        FREQDOMAINVALUES[id]['Hdu'] = np.divide(FREQDOMAINVALUES[id]['INPS'][IX],FREQDOMAINVALUES[id]['DISTS'][IX])
    return timedomainvalues, FREQDOMAINVALUES

def plotThings(M,Rid,Did,Bothid,Evenid,Oddid,IX,IXE,IXO,subject):
    # find ffts
    ref, REF = findFFT(Rid,IX,subject)
    dis,DIS = findFFT(Did,IX,subject)
    both,BOTH = findFFT(Bothid,IX,subject)


#   UNCOMMENT to plot Hdu, Hru for all trials
#     plt.figure(figsize=(20,10))
#     # plot Hru
#     plt.subplot(241)
#     linecolor = 0.2
#     for id in Rid:
#         plt.loglog(freq[IX],np.absolute(REF[id]['Hru']),color='%f'%linecolor)
#         linecolor += .08
#     plt.ylabel('Magnitude')
#     plt.title('Hru - ref only')

#     plt.subplot(245)
#     linecolor = 0.2
#     for id in Rid:
#         plt.semilogx(freq[IX],np.angle(REF[id]['Hru'])*180/m.pi,color='%f'%linecolor)
#         linecolor += .08
#     plt.ylabel('Phase - deg')
#     plt.xlabel('Frequency - Hz')

#     # plot Hdu
#     plt.subplot(242)
#     linecolor = 0.2
#     for id in Did:
#         plt.loglog(freq[IX],np.absolute(DIS[id]['Hdu']),color='%f'%linecolor)
#         linecolor += .08
#     plt.ylabel('Magnitude')
#     plt.title('Hdu - dis only')

#     plt.subplot(246)
#     linecolor = 0.2
#     for id in Did:
#         plt.semilogx(freq[IX],np.angle(DIS[id]['Hdu'])*180/m.pi,color='%f'%linecolor)
#         linecolor += .08
#     plt.ylabel('Phase - deg')
#     plt.xlabel('Frequency - Hz')

#     # plot BOTH
#     plt.subplot(243)
#     linecolor=.2
#     for id in Evenid:
#         plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hru'][0::2]),color='%f'%linecolor)
#         linecolor+=.08
#     for id in Oddid:
#         plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hru'][1::2]))
#     plt.ylabel('Magnitude')
#     plt.title('Hru - BOTH')

#     plt.subplot(247)
#     linecolor=.2
#     for id in Evenid:
#         plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hru'][0::2])*180/m.pi,color='%f'%linecolor)
#         linecolor+=.08
#     for id in Oddid:
#         plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hru'][1::2])*180/m.pi)
#     plt.ylabel('Phase - deg')
#     plt.xlabel('Frequency - Hz')

#     # plot Hdu
#     plt.subplot(244)
#     linecolor=.2
#     for id in Evenid:
#         plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hdu'][1::2]))
#     for id in Oddid:
#         plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hdu'][0::2]),color='%f'%linecolor)
#         linecolor+=.08
#     plt.ylabel('Magnitude')
#     plt.title('Hdu - both')

#     plt.subplot(248)
#     linecolor=.2
#     for id in Evenid:
#         plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hdu'][1::2])*180/m.pi)
#     for id in Oddid:
#         plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hdu'][0::2])*180/m.pi,color='%f'%linecolor)
#         linecolor+=.08
#     plt.ylabel('Phase - deg')
#     plt.xlabel('Frequency - Hz')



    # now find averages for Hru and Hdu for R AND D SOLO
    HruRavg = np.zeros((len(IX),))
    for id in Rid:
        HruRavg = HruRavg+REF[id].get('Hru')
    HruRavg = HruRavg/len(REF)

    HduDavg = np.zeros((len(IX),))
    for id in Did:
        HduDavg = HduDavg+DIS[id].get('Hdu')
    HduDavg = HduDavg/len(DIS)

    # and get averages for FOB and FOF
    BD = np.divide(-HduDavg,(np.multiply(M,(1.0+0j)+HduDavg)))
    FR = np.multiply(HruRavg,(1.0+0j)+np.multiply(BD,M))-BD

    BonlyR = np.divide(HruRavg,(1.0+0j)-np.multiply(HruRavg,M))
    BonlyD = np.divide(-HduDavg,(np.multiply(M,(1.0+0j)+HduDavg)))

    # UNCOMMENT to get individual plots for feedback and feedforward policies from R only and D only trials
    # plt.figure()
    # plt.subplot(221)
    # num = len(Rid)
    # plt.title('Estimated Feedback/Feedforward Policies from R, D trials (n=%i)'%num)
    # plt.loglog(freq[IX],np.absolute(BD),'b',label='Feedback')
    # plt.loglog(freq[IX],np.absolute(FR),'r',label='Feedforward')
    # plt.loglog(freq[IX],np.absolute(np.divide(1,M)),'r--',label='Model Inverse Prediction')
    # plt.legend()
    # plt.ylabel('Gain')
    #
    # plt.subplot(223)
    # plt.semilogx(freq[IX],np.angle(BD)*180/m.pi,'b',label='Feedback')
    # plt.semilogx(freq[IX],np.angle(FR)*180/m.pi,'r',label='Feedforward')
    # plt.semilogx(freq[IX],np.angle(np.divide(1,M))*180/m.pi,'r--',label='Model Inverse Prediction')
    # plt.legend()
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Freq - Hz')

    # get averages for BOTH
    # this is for Hru
    tempE = np.zeros((len(IXE),))
    for id in Evenid:
        tempE = tempE +BOTH[id].get('Hru')[0::2]
    tempE = tempE/len(Evenid)

    tempO = np.zeros((len(IXO),))
    for id in Oddid:
        tempO = tempO +BOTH[id].get('Hru')[1::2]
    tempO = tempO/len(Oddid)

    HruRDavg = np.zeros((len(IX),),dtype=complex)
    HruRDavg[::2] = tempE
    HruRDavg[1::2] = tempO

    # this is for Hdu
    tempO = np.zeros((len(IXO),))
    for id in Evenid:
        tempO = tempO +BOTH[id].get('Hdu')[1::2]
    tempO = tempO/len(Evenid)

    tempE = np.zeros((len(IXE),))
    for id in Oddid:
        tempE = tempE +BOTH[id].get('Hdu')[::2]
    tempE = tempE/len(Oddid)

    HduRDavg = np.zeros((len(IX),),dtype=complex)
    HduRDavg[1::2] = tempO
    HduRDavg[::2] = tempE

    # find F and B
    BRD = np.divide(-HduRDavg,(np.multiply(M,(1.0+0j)+HduRDavg)))
    FRD = np.multiply(HruRDavg,(1.0+0j)+np.multiply(BRD,M))-BRD

    # UNCOMMENT to get individual plots for feedback and feedforward policies for R+D trials
    # plt.subplot(222)
    # num = len(Oddid)
    # plt.title('Estimated Feedback/Feedforward Policies from R+D trials (n=%i)'%num)
    # plt.loglog(freq[IX],np.absolute(BRD),'b',label='Feedback')
    # plt.loglog(freq[IX],np.absolute(FRD),'r',label='Feedforward')
    # plt.loglog(freq[IX],np.absolute(np.divide(1.0,M)),'r--',label='Model Inverse Prediction')
    # plt.legend()
    # plt.ylabel('Gain')
    #
    # plt.subplot(224)
    # plt.semilogx(freq[IX],np.angle(BRD)*180/m.pi,'b',label='Feedback')
    # plt.semilogx(freq[IX],np.angle(FRD)*180/m.pi,'r',label='Feedforward')
    # plt.semilogx(freq[IX],np.angle(np.divide(1.0,M))*180/m.pi,'r--',label='Model Inverse Prediction')
    # plt.legend()
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Freq - Hz')
    return BD,FR,BRD,FRD,dis,DIS,ref,REF,both,BOTH,HruRavg,HduDavg,HruRDavg,HduRDavg,BonlyD,BonlyR
