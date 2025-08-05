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
def findFFT(listofIDs,IX):
    trials = dict()
    for id in listofIDs:
      fis = sorted(glob.glob(os.path.join('data',subject,id+'.npz')))#'*_'+protocol+'_'+id+'.npz')))
      if len(fis) > 1:
        dbg('WARNING -- repeated trials for id ='+id)
      assert len(fis) > 0, 'ERROR -- no data for id ='+id
      fi = fis[-1]
      print('LOAD '+fi)
      trial = dict(np.load(fi,encoding="latin1"))
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
        FREQDOMAINVALUES[id]['OUTS'] = np.fft.fft(timedomainvalues[id]['outs'])[:int(timedomainvalues[id]['outs'].shape[0]/2)]/TEMPFREQ
        FREQDOMAINVALUES[id]['INPS'] = np.fft.fft(timedomainvalues[id]['inps'])[:int(timedomainvalues[id]['inps'].shape[0]/2)]/TEMPFREQ
        FREQDOMAINVALUES[id]['REFS'] = np.fft.fft(timedomainvalues[id]['refs'])[:int(timedomainvalues[id]['refs'].shape[0]/2)]/TEMPFREQ
        FREQDOMAINVALUES[id]['DISTS'] = np.fft.fft(timedomainvalues[id]['dists'])[:int(timedomainvalues[id]['dists'].shape[0]/2)]/TEMPFREQ
        FREQDOMAINVALUES[id]['Hru'] = np.divide(FREQDOMAINVALUES[id]['INPS'][IX],FREQDOMAINVALUES[id]['REFS'][IX])
        FREQDOMAINVALUES[id]['Hdu'] = np.divide(FREQDOMAINVALUES[id]['INPS'][IX],FREQDOMAINVALUES[id]['DISTS'][IX])
    return timedomainvalues, FREQDOMAINVALUES

def plotThings(M,Rid,Did,Bothid,Evenid,Oddid,IX,IXE,IXO):
    # find ffts
    ref, REF = findFFT(Rid,IX)
    dis,DIS = findFFT(Did,IX)
    both,BOTH = findFFT(Bothid,IX)

    # making figures
    # plot Hdu, Hru for all trials
    # plt.figure()
    # # plot Hru
    # plt.subplot(241)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.loglog(freq[IX],np.absolute(REF[id]['Hru']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Magnitude')
    # plt.title('Hru - ref only')
    #
    # plt.subplot(245)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.semilogx(freq[IX],np.angle(REF[id]['Hru'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot Hdu
    # plt.subplot(242)
    # linecolor = 0.2
    # for id in Did:
    #     plt.loglog(freq[IX],np.absolute(DIS[id]['Hdu']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Magnitude')
    # plt.title('Hdu - dis only')
    #
    # plt.subplot(246)
    # linecolor = 0.2
    # for id in Did:
    #     plt.semilogx(freq[IX],np.angle(DIS[id]['Hdu'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot BOTH
    # plt.subplot(243)
    # linecolor=.2
    # for id in Evenid:
    #     plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hru'][0::2]),color='%f'%linecolor)
    #     linecolor+=.08
    # for id in Oddid:
    #     plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hru'][1::2]))
    # plt.ylabel('Magnitude')
    # plt.title('Hru - BOTH')
    #
    # plt.subplot(247)
    # linecolor=.2
    # for id in Evenid:
    #     plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hru'][0::2])*180/m.pi,color='%f'%linecolor)
    #     linecolor+=.08
    # for id in Oddid:
    #     plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hru'][1::2])*180/m.pi)
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot Hdu
    # plt.subplot(244)
    # linecolor=.2
    # for id in Evenid:
    #     plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hdu'][1::2]))
    # for id in Oddid:
    #     plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hdu'][0::2]),color='%f'%linecolor)
    #     linecolor+=.08
    # plt.ylabel('Magnitude')
    # plt.title('Hdu - both')
    #
    # plt.subplot(248)
    # linecolor=.2
    # for id in Evenid:
    #     plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hdu'][1::2])*180/m.pi)
    # for id in Oddid:
    #     plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hdu'][0::2])*180/m.pi,color='%f'%linecolor)
    #     linecolor+=.08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')


    # now find averages for Hru and Hdu for R AND D SOLO


    HruRall = np.zeros((len(Rid),len(IX)),dtype=np.complex_)
    for num,id in enumerate(Rid):
        HruRall[num,:] = REF[id].get('Hru')
    HruR25 = {}
    HruR75 = {}
    HruRavg = {}
    HruR25['mag'],HruRavg['mag'],HruR75['mag'] = np.percentile(np.absolute(HruRall), [25,50,75],axis=0)
    HruR25['ang'],HruRavg['ang'],HruR75['ang'] = np.percentile(abs(np.angle(HruRall)), [25,50,75],axis=0)

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

    # plt.figure()
    plt.figure()
    plt.subplot(121)
    plt.fill_between(freq[IX], HruR25['mag'], HruR75['mag'], color=(0.1,0.1,0.4))#"#3F5D7D")
    plt.semilogx(freq[IX],HruRavg['mag'],color=(.6,.6,1),label='r all')
    plt.loglog(freq[IXE],np.absolute(tempE),'bo',label='r Even')
    plt.loglog(freq[IXO],np.absolute(tempO),'go',label='r Odd')
    plt.xlabel('freq (Hz)')
    plt.ylabel('Magnitude')

    plt.subplot(122)
    plt.fill_between(freq[IX], HruR25['ang'],HruR75['ang'], color=(0.4,0.1,0.1))#"#3F5D7D")
    plt.semilogx(freq[IX],HruRavg['ang'],color=(1,.6,.6),label='r all')
    plt.semilogx(freq[IXE],abs(np.angle(tempE)),'bo',label='r Even')
    plt.semilogx(freq[IXO],abs(np.angle(tempO)),'go',label='r Odd')
    plt.xlabel('freq (Hz)')
    plt.ylabel('angle (rad)')
    plt.legend()

    # for id in Rid:
    #     HruRavg = HruRavg+REF[id].get('Hru')
    # HruRavg = HruRavg/len(REF)

    # # plot Hru
    # plt.figure()
    # plt.subplot(121)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.loglog(freq[IX],np.absolute(REF[id]['Hru']),color='%f'%linecolor,marker='o')
    #     linecolor += .08
    # plt.xlabel('freq (Hz)')
    # plt.ylabel('Magnitude')
    # plt.title('Hru - ref only')
    #
    # plt.subplot(122)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.semilogx(freq[IX],np.angle(REF[id]['Hru'])*180/m.pi,color='%f'%linecolor,marker='o')
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # HduDavg = np.zeros((len(IX),))
    # for id in Did:
    #     HduDavg = HduDavg+DIS[id].get('Hdu')
    # HduDavg = HduDavg/len(DIS)
    # # and get averages for FOB and FOF
    # BD = np.divide(-HduDavg,(np.multiply(M,(1.0+0j)+HduDavg)))
    # FR = np.multiply(HruRavg,(1.0+0j)+np.multiply(BD,M))-BD
    #
    # for id in Rid:
    #     REF[id]['FR'] = np.multiply(REF[id].get('Hru'),(1.0+0j)+np.multiply(BD,M))-BD

    # # plot Hru
    # plt.figure()
    # plt.subplot(121)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.loglog(freq[IX],np.absolute(REF[id]['FR']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.loglog(freq[IX],np.absolute(np.divide(1.0,M)),'r--',label='Model Inverse Prediction')
    # plt.ylabel('Magnitude')
    # plt.title('FF - ref only')
    #
    # plt.subplot(122)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.semilogx(freq[IX],np.angle(REF[id]['FR'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.semilogx(freq[IX],np.angle(np.divide(1.0,M))*180/m.pi,'r--',label='Model Inverse Prediction')
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
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

    # # get averages for BOTH
    # # this is for Hdu
    # tempO = np.zeros((len(IXO),))
    # for id in Evenid:
    #     tempO = tempO +BOTH[id].get('Hdu')[1::2]
    # tempO = tempO/len(Evenid)
    #
    # tempE = np.zeros((len(IXE),))
    # for id in Oddid:
    #     tempE = tempE +BOTH[id].get('Hdu')[::2]
    # tempE = tempE/len(Oddid)
    #
    # HduRDavg = np.zeros((len(IX),),dtype=complex)
    # HduRDavg[1::2] = tempO
    # HduRDavg[::2] = tempE
    # BDBOTH = np.divide(-HduRDavg,(np.multiply(M,(1.0+0j)+HduRDavg)))

    # get Hru
    # tempE = np.zeros((len(IXE),))
    # for id in Evenid:
    #     tempE = tempE +BOTH[id].get('Hru')[0::2]
    # tempE = tempE/len(Evenid)
    #
    # tempO = np.zeros((len(IXO),))
    # for id in Oddid:
    #     tempO = tempO +BOTH[id].get('Hru')[1::2]
    # tempO = tempO/len(Oddid)
    #
    # HruRDavg = np.zeros((len(IX),),dtype=complex)
    # HruRDavg[::2] = tempE
    # HruRDavg[1::2] = tempO
    #
    # plt.figure()
    # plt.title('f1 and f2')
    # plt.subplot(121)
    # plt.loglog(freq[IXE],np.absolute(tempE),'bo',label='r Even')
    # plt.loglog(freq[IXO],np.absolute(tempO),'go',label='r Odd')
    # plt.xlabel('freq (Hz)')
    # plt.ylabel('Magnitude')
    #
    # plt.subplot(122)
    # plt.semilogx(freq[IXE],np.angle(tempE),'bo',label='r Even')
    # plt.semilogx(freq[IXO],np.angle(tempO),'go',label='r Odd')
    # plt.xlabel('freq (Hz)')
    # plt.ylabel('angle (rad)')
    # plt.legend()
    #

    #
    # tempE = {}#np.zeros((len(IXE),))
    # for id in Evenid:
    #     tempE[id] = np.multiply(BOTH[id].get('Hru')[0::2],(1.0+0j)+np.multiply(BDBOTH[0::2],M[0::2]))-BDBOTH[0::2]
    #
    # tempO = {}
    # for id in Oddid:
    #     tempO[id] = np.multiply(BOTH[id].get('Hru')[1::2],(1.0+0j)+np.multiply(BDBOTH[1::2],M[1::2]))-BDBOTH[1::2]
    #
    # FRBOTH = np.zeros((min(len(Evenid),len(Oddid)),len(M)),dtype=np.complex_)
    # for num in range(0,min(len(Evenid),len(Oddid))):
    #     FRBOTH[num,::2] = tempE[Evenid[num]]
    #     FRBOTH[num,1::2]= tempO[Oddid[num]]

    # plt.figure()
    # plt.subplot(121)
    # linecolor = 0.2
    # for num in range(0,len(FRBOTH)):
    #     plt.loglog(freq[IX],np.absolute(FRBOTH[num]),color='%f'%linecolor)
    #     linecolor += .08
    # plt.loglog(freq[IX],np.absolute(np.divide(1.0,M)),'r--',label='Model Inverse Prediction')
    # plt.ylabel('Magnitude')
    # plt.title('FF - both presented')
    #
    # plt.subplot(122)
    # linecolor = 0.2
    # for num in range(0,len(FRBOTH)):
    #     plt.semilogx(freq[IX],np.angle(FRBOTH[num])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.semilogx(freq[IX],np.angle(np.divide(1.0,M))*180/m.pi,'r--',label='Model Inverse Prediction')
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')

    # HruRDavg = np.zeros((len(IX),),dtype=complex)
    # HruRDavg[::2] = tempE
    # HruRDavg[1::2] = tempO
    #
    #
    # # find F and B
    # BRD = np.divide(-HduRDavg,(np.multiply(M,(1.0+0j)+HduRDavg)))
    # FRD = np.multiply(HruRDavg,(1.0+0j)+np.multiply(BRD,M))-BRD
    #
    # for id in Rid:
    #     REF[id]['FR'] = np.multiply(REF[id].get('Hru'),(1.0+0j)+np.multiply(BD,M))-BD
    #
    # # plot Hru
    # plt.figure()
    # plt.subplot(121)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.loglog(freq[IX],np.absolute(REF[id]['FR']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Magnitude')
    # plt.title('FF - ref only')
    #
    # plt.subplot(122)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.semilogx(freq[IX],np.angle(REF[id]['FR'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
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
    #
    #
    # plt.figure()
    # plt.plot(ref[Rid[0]]['times'],ref[Rid[0]]['refs'],label='Reference')
    # plt.plot(ref[Rid[0]]['times'],ref[Rid[0]]['outs'],label='Output')
    # plt.xlabel('time')
    # plt.ylabel('distance from 0')
    # plt.legend()



## TODO CHANGE TO PLOTTING EACH FF AND FB
    # plt.figure()
    # # plot Hru
    # plt.subplot(241)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.loglog(freq[IX],np.absolute(REF[id]['Hru']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Magnitude')
    # plt.title('Hru - ref only')
    #
    # plt.subplot(245)
    # linecolor = 0.2
    # for id in Rid:
    #     plt.semilogx(freq[IX],np.angle(REF[id]['Hru'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot Hdu
    # plt.subplot(242)
    # linecolor = 0.2
    # for id in Did:
    #     plt.loglog(freq[IX],np.absolute(DIS[id]['Hdu']),color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Magnitude')
    # plt.title('Hdu - dis only')
    #
    # plt.subplot(246)
    # linecolor = 0.2
    # for id in Did:
    #     plt.semilogx(freq[IX],np.angle(DIS[id]['Hdu'])*180/m.pi,color='%f'%linecolor)
    #     linecolor += .08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot BOTH
    # plt.subplot(243)
    # linecolor=.2
    # for id in Evenid:
    #     plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hru'][0::2]),color='%f'%linecolor)
    #     linecolor+=.08
    # for id in Oddid:
    #     plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hru'][1::2]))
    # plt.ylabel('Magnitude')
    # plt.title('Hru - BOTH')
    #
    # plt.subplot(247)
    # linecolor=.2
    # for id in Evenid:
    #     plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hru'][0::2])*180/m.pi,color='%f'%linecolor)
    #     linecolor+=.08
    # for id in Oddid:
    #     plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hru'][1::2])*180/m.pi)
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # # plot Hdu
    # plt.subplot(244)
    # linecolor=.2
    # for id in Evenid:
    #     plt.loglog(freq[IXO],np.absolute(BOTH[id]['Hdu'][1::2]))
    # for id in Oddid:
    #     plt.loglog(freq[IXE],np.absolute(BOTH[id]['Hdu'][0::2]),color='%f'%linecolor)
    #     linecolor+=.08
    # plt.ylabel('Magnitude')
    # plt.title('Hdu - both')
    #
    # plt.subplot(248)
    # linecolor=.2
    # for id in Evenid:
    #     plt.semilogx(freq[IXO],np.angle(BOTH[id]['Hdu'][1::2])*180/m.pi)
    # for id in Oddid:
    #     plt.semilogx(freq[IXE],np.angle(BOTH[id]['Hdu'][0::2])*180/m.pi,color='%f'%linecolor)
    #     linecolor+=.08
    # plt.ylabel('Phase - deg')
    # plt.xlabel('Frequency - Hz')
    #
    # plt.figure()
    # plt.scatter(FR.real,FR.imag,label='just R')
    # plt.scatter(FRD.real,FRD.imag, label='R+D')
    # plt.legend()

    #plt.scatter(REF[Rid[-1]]['Hru'].real,REF[Rid[-1]]['Hru'].imag)

help = """
usage:
  analysis subject protocol [formats]

add filename formats to save files in those formats:
  visualization subject protocol png,mp4
"""

def dbg(s):
  print (s)
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


# importing protocol (su17v3)
#proto = importlib.import_module('protocols.'+protocol)

# find files for specific subject
#fis = glob.glob(os.path.join('data',subject,'*_'+protocol+'*'))
#ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
#ids = [id for id in ids if ('.csv') not in id[-4:]]
sys.path.insert(0,'C:\\Users\\BRL\\Momona\'s Google Drive\\Yamagami Lab_\\Slider Project\\hcps\\develop\\Analysis\\Momona\'s code\\protocols')
import su17v3 as proto
#proto = importlib.import_module('protocols.'+protocol)
sys.path.insert(0,'C:\\Users\\BRL\\Momona\'s Google Drive\\Yamagami Lab_\\Slider Project\\hcps\\develop')
# find files for specific subject
fis = glob.glob('C:/Users/BRL/Momona\'s Google Drive/Yamagami Lab_/Slider Project/hcps/develop/Analysis/Momona\'s code/data/'+subject+'/*.npz')#.npz')
fis = [os.path.basename(fi) for fi in fis]
#fis = glob.glob(os.path.join('data'+subject,'*_'+subject+'*'))
#ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
ids = sorted(list(set([fi.strip('.npz') for fi in fis])))
ids = [id for id in ids if ('st2') not in id]
ids = [id for id in ids if ('st1') not in id]
foids = [id for id in ids if ('fo') in id]
soids = [id for id in ids if ('fo') not in id]
foRid = [id for id in foids if ('d-zer') in id]
soRid = [id for id in soids if ('d-zer') in id]
foDid = [id for id in foids if ('r-zer') in id]
soDid = [id for id in soids if ('r-zer') in id]
foBothid = [id for id in foids if ('zer') not in id]
foEvenid = [id for id in foBothid if ('r-sos+E') in id or ('r-sos-E') in id]
foOddid = [id for id in foBothid if ('r-sos-O') in id or ('r-sos+O') in id]
soBothid = [id for id in soids if ('zer') not in id]
soEvenid = [id for id in soBothid if ('r-sos+E') in id or ('r-sos-E') in id]
soOddid = [id for id in soBothid if ('r-sos-O') in id or ('r-sos+O') in id]

# Define
primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31])# max =37
IX = primes*2
primesE = np.asarray([2,5,11,17,23,31])#primes[::2]
primesO = np.asarray([3,7,13,19,29])#primes[1::2]
IXE = primesE*2
IXO = primesO*2
freq = np.linspace(0,599.5,1200)#refs.shape[0]/2)
tempfreq = np.ones((2400,))
TEMPFREQ = np.fft.fft(tempfreq)[0]
M = 1/(1j*(2*m.pi*freq[IX]))*20
M = np.asarray(M)
# first order
plotThings(M,foRid,foDid,foBothid,foEvenid,foOddid,IX,IXE,IXO)

# # second order
# # Define
# primes = np.asarray([2, 3, 5, 7, 11, 13,17])# max =37
# IX = primes*2
# primesE = np.asarray([2,5,11,17])#primes[::2]
# primesO = np.asarray([3,7,13])#primes[1::2]
# IXE = primesE*2
# IXO = primesO*2
# freq = np.linspace(0,599.5,1200)#refs.shape[0]/2)
# tempfreq = np.ones((2400,))
# TEMPFREQ = np.fft.fft(tempfreq)[0]
# s = 1.0j*2*m.pi*freq[IX]/20
# M = 1.0/(np.multiply(s,s)+s)#(1.0/(-(2*m.pi*freq[IX])**2+1j*(2*m.pi*freq[IX])))*20
# M = np.asarray(M)
# plotThings(M,soRid,soDid,soBothid,soEvenid,soOddid,IX,IXE,IXO)

plt.show()
