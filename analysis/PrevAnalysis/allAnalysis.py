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
plt.close()

    # Define
primes = np.asarray([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31])# max =37
foIX = primes*2
primesE = np.asarray([2,5,11,17,23,31])#primes[::2]
primesO = np.asarray([3,7,13,19,29])#primes[1::2]
foIXE = primesE*2
foIXO = primesO*2
freq = np.linspace(0,599.5,1200)#refs.shape[0]/2)
tempfreq = np.ones((2400,))
TEMPFREQ = np.fft.fft(tempfreq)[0]
foM = 1/(1j*(2*m.pi*freq[foIX]))*20
foM = np.asarray(foM)

# Define
primes = np.asarray([2, 3, 5, 7, 11, 13,17])# max =37
soIX = primes*2
primesE = np.asarray([2,5,11,17])#primes[::2]
primesO = np.asarray([3,7,13])#primes[1::2]
soIXE = primesE*2
soIXO = primesO*2
s = 1.0j*2*m.pi*freq[soIX]/20
soM = 1.0/(np.multiply(s,s)+s)#(1.0/(-(2*m.pi*freq[IX])**2+1j*(2*m.pi*freq[IX])))*20
soM = np.asarray(soM)

#
def findFFT(listofIDs,IX):
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
    return BD,FR,BRD,FRD,dis,DIS,ref,REF,both,BOTH



help = """
usage:
  analysis subject protocol [formats]

add filename formats to save files in those formats:
  visualization subject protocol png,mp4
"""

def dbg(s):
  print s
# CHECK THAT ALL RELEVANT PARAMETERS ARE SPECIFIED IN COMMAND
args = sys.argv
# if len(args) < 2:
#   dbg('\nABORT -- no subject specified')
#   dbg(help)
#   sys.exit(0)
# if len(args) < 3:
#   dbg('\nABORT -- no protocol specified')
#   dbg(help)
#   sys.exit(0)
# subject = args[1]
# protocol = args[2]
# fmts = []
# if len(args) >= 4:
#   fmts = args[3].split(',')
#

# importing protocol (su17v3)
protocol = 'su17v3'
proto = importlib.import_module('protocols.'+protocol)

def takeAverage(subjects):
    # find files for specific subject
    global subject
    subject = subjects
    fis = glob.glob(os.path.join('data',subject,'*_'+protocol+'*'))
    ids = sorted(list(set([fi.strip('.npz').split('_',2)[2] for fi in fis])))
    ids = [id for id in ids if ('.csv') not in id[-4:]]
    ids = [id for id in ids if id[-3:] not in ['rej','oob','rst','max','st1','st2']]
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


    # first order
    BD,FR,BRD,FRD,dis,DIS,ref,REF,both,BOTH = plotThings(foM,foRid,foDid,foBothid,foEvenid,foOddid,foIX,foIXE,foIXO)

    # second order

    #BD,FR,BRD,FRD = plotThings(soM,soRid,soDid,soBothid,soEvenid,soOddid,soIX,soIXE,soIXO)
    return BD,FR,BRD,FRD,dis,DIS,ref,REF,both,BOTH

subjects = np.asarray(['ef','dp','gy','iz','ly','wv','zu'])
print subjects.shape

data = {}
dataAvg = {}
dataAvg['BDavg'] = np.zeros((11,len(subjects)))
dataAvg['FRavg'] = np.zeros((11,len(subjects)))
dataAvg['BRDavg'] =np.zeros((11,len(subjects)))
dataAvg['FRDavg'] =np.zeros((11,len(subjects)))

angleAvg = {}
angleAvg['BDavg'] = np.zeros((11,len(subjects)))
angleAvg['FRavg'] =np.zeros((11,len(subjects)))
angleAvg['BRDavg'] =np.zeros((11,len(subjects)))
angleAvg['FRDavg'] = np.zeros((11,len(subjects)))
numcount = 0;
for sub in subjects:
    data[sub] = {}
    data[sub]['BD'],data[sub]['FR'],data[sub]['BRD'],data[sub]['FRD'],
        data[sub]['dis'],data[sub]['DIS'],data[sub]['ref'],data[sub]['REF'],
        data[sub]['both'],data[sub]['BOTH']= takeAverage(sub)
    data[sub]['BDabs'] = np.absolute(data[sub]['BD'])
    data[sub]['FRabs'] = np.absolute(data[sub]['FR'])
    data[sub]['BRDabs'] = np.absolute(data[sub]['BRD'])
    data[sub]['FRDabs'] = np.absolute(data[sub]['FRD'])

    data[sub]['BDangle'] = np.angle(data[sub]['BD'])*180/m.pi
    data[sub]['FRangle'] = np.angle(data[sub]['FR'])*180/m.pi
    data[sub]['BRDangle'] = np.angle(data[sub]['BRD'])*180/m.pi
    data[sub]['FRDangle'] = np.angle(data[sub]['FRD'])*180/m.pi

    dataAvg['BDavg'][:,numcount]=data[sub]['BDabs']
    dataAvg['FRavg'][:,numcount]=data[sub]['FRabs']
    dataAvg['BRDavg'][:,numcount]=data[sub]['BRDabs']
    dataAvg['FRDavg'][:,numcount]=data[sub]['FRDabs']

    angleAvg['BDavg'][:,numcount]=np.angle(data[sub]['BD'])*180/m.pi
    angleAvg['FRavg'][:,numcount]=np.angle(data[sub]['FR'])*180/m.pi
    angleAvg['BRDavg'][:,numcount]=np.angle(data[sub]['BRD'])*180/m.pi
    angleAvg['FRDavg'][:,numcount]=np.angle(data[sub]['FRD'])*180/m.pi
    numcount +=1

#find stdev and mean for magnitude
dataAvg['BDallavg'] = np.mean(dataAvg['BDavg'],1)
dataAvg['FRallavg'] = np.mean(dataAvg['FRavg'],1)
dataAvg['BRDallavg'] = np.mean(dataAvg['BRDavg'],1)
dataAvg['FRDallavg'] = np.mean(dataAvg['FRDavg'],1)

dataAvg['BDallstd'] = np.std(dataAvg['BDavg'],1)
dataAvg['FRallstd'] = np.std(dataAvg['FRavg'],1)
dataAvg['BRDallstd'] = np.std(dataAvg['BRDavg'],1)
dataAvg['FRDallstd'] = np.std(dataAvg['FRDavg'],1)

# find mean and std for angle
angleAvg['BDallavg'] = np.mean(angleAvg['BDavg'],1)
angleAvg['FRallavg'] = np.mean(angleAvg['FRavg'],1)
angleAvg['BRDallavg'] = np.mean(angleAvg['BRDavg'],1)
angleAvg['FRDallavg'] = np.mean(angleAvg['FRDavg'],1)

angleAvg['BDallstd'] = np.std(angleAvg['BDavg'],1)
angleAvg['FRallstd'] = np.std(angleAvg['FRavg'],1)
angleAvg['BRDallstd'] = np.std(angleAvg['BRDavg'],1)
angleAvg['FRDallstd'] = np.std(angleAvg['FRDavg'],1)

# dataAvg['BDavg'] = dataAvg['BDavg']/len(subjects)
# dataAvg['FRavg'] = dataAvg['FRavg']/len(subjects)
# dataAvg['BRDavg'] = dataAvg['BRDavg']/len(subjects)
# dataAvg['FRDavg'] = dataAvg['FRDavg']/len(subjects)
#
# angleAvg['BDavg'] = angleAvg['BDavg']/len(subjects)
# angleAvg['FRavg'] = angleAvg['FRavg']/len(subjects)
# angleAvg['BRDavg'] = angleAvg['BRDavg']/len(subjects)
# angleAvg['FRDavg'] = angleAvg['FRDavg']/len(subjects)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

plt.figure()
plt.subplot(221)
#plt.title('Estimated Feedback/Feedforward Policies from R, D trials')
plt.fill_between(freq[foIX], dataAvg['BDallavg'] - dataAvg['BDallstd'],
                 dataAvg['BDallavg'] + dataAvg['BDallstd'], color=(0.1,0.1,0.4))#"#3F5D7D")
plt.loglog(freq[foIX],dataAvg['BDallavg'],color=(.6,.6,1),label='Feedback')

plt.fill_between(freq[foIX], dataAvg['FRallavg'] - dataAvg['FRallstd'],
                 dataAvg['FRallavg'] + dataAvg['FRallstd'], color=(0.4,0.1,0.1))#"#3F5D7D")
plt.loglog(freq[foIX],dataAvg['FRallavg'],color=(1,.6,.6),label='Feedforward')
plt.loglog(freq[foIX],np.absolute(np.divide(1,foM)),color=(1,.6,.6),linestyle='--',label='Model Inverse Prediction')
#plt.legend()
plt.ylabel('Gain')
#plt.xlabel('Frequency (Hz)')


plt.subplot(223)
#plt.title('Estimated Feedback/Feedforward Policies from R, D trials')
plt.fill_between(freq[foIX], angleAvg['BDallavg'] - angleAvg['BDallstd'],
                 angleAvg['BDallavg'] + angleAvg['BDallstd'], color=(0.1,0.1,0.4))#"#3F5D7D")
plt.semilogx(freq[foIX],angleAvg['BDallavg'],color=(.6,.6,1),label='Feedback')

plt.fill_between(freq[foIX], angleAvg['FRallavg'] - angleAvg['FRallstd'],
                 angleAvg['FRallavg'] + angleAvg['FRallstd'], color=(0.4,0.1,0.1))#"#3F5D7D")
plt.semilogx(freq[foIX],angleAvg['FRallavg'],color=(1,.6,.6),label='Feedforward')
plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,color=(1,.6,.6),linestyle='--',label='Model Inverse Prediction')
#plt.legend()
plt.ylabel('Angle')
plt.xlabel('Frequency (Hz)')


plt.subplot(222)
#plt.title('Estimated Feedback/Feedforward Policies from R+D (interleaved) trials')
plt.fill_between(freq[foIX], dataAvg['BRDallavg'] - dataAvg['BRDallstd'],
                 dataAvg['BRDallavg'] + dataAvg['BRDallstd'], color=(0.1,0.1,0.4))#"#3F5D7D")
plt.loglog(freq[foIX],dataAvg['BRDallavg'],color=(.6,.6,1),label='Feedback')

plt.fill_between(freq[foIX], dataAvg['FRDallavg'] - dataAvg['FRDallstd'],
                 dataAvg['FRDallavg'] + dataAvg['FRDallstd'], color=(0.4,0.1,0.1))#"#3F5D7D")
plt.loglog(freq[foIX],dataAvg['FRDallavg'],color=(1,.6,.6),label='Feedforward')
plt.loglog(freq[foIX],np.absolute(np.divide(1,foM)),color=(1,.6,.6),linestyle='--',label='Model Inverse Prediction')
#plt.legend()
plt.ylabel('Gain')
#plt.xlabel('Frequency (Hz)')

plt.subplot(224)
#plt.title('Estimated Feedback/Feedforward Policies from R+D (interleaved) trials')
plt.fill_between(freq[foIX], angleAvg['BRDallavg'] - angleAvg['BRDallstd'],
                 angleAvg['BRDallavg'] + angleAvg['BRDallstd'], color=(0.1,0.1,0.4))#"#3F5D7D")
plt.semilogx(freq[foIX],angleAvg['BRDallavg'],color=(.6,.6,1),label='Feedback')

plt.fill_between(freq[foIX], angleAvg['FRDallavg'] - angleAvg['FRDallstd'],
                 angleAvg['FRDallavg'] + angleAvg['FRDallstd'], color=(0.4,0.1,0.1))#"#3F5D7D")
plt.semilogx(freq[foIX],angleAvg['FRDallavg'],color=(1,.6,.6),label='Feedforward')
plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,color=(1,.6,.6),linestyle='--',label='Model Inverse Prediction')
#plt.legend()
plt.ylabel('Angle')
#
# plt.subplot(223)
# plt.semilogx(freq[foIX],angleAvg['BDavg'],'b',label='Feedback')
# plt.semilogx(freq[foIX],angleAvg['FRavg'],'r',label='Feedforward')
# plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,'r--',label='Model Inverse Prediction')
# plt.legend()
# plt.ylabel('Phase - deg')
# plt.xlabel('Freq - Hz')
#
#
# plt.subplot(222)
# plt.title('Estimated Feedback/Feedforward Policies from R+D trials')
# plt.loglog(freq[foIX],dataAvg['BRDavg'],'b',label='Feedback')
# plt.loglog(freq[foIX],dataAvg['FRDavg'],'r',label='Feedforward')
# plt.loglog(freq[foIX],np.absolute(np.divide(1,foM)),'r--',label='Model Inverse Prediction')
# plt.legend()
# plt.ylabel('Gain')
#
# plt.subplot(224)
# plt.semilogx(freq[foIX],angleAvg['BRDavg'],'b',label='Feedback')
# plt.semilogx(freq[foIX],angleAvg['FRDavg'],'r',label='Feedforward')
# plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,'r--',label='Model Inverse Prediction')
# plt.legend()
# plt.ylabel('Phase - deg')
# plt.xlabel('Freq - Hz')



plt.figure()
plt.subplot(221)
#plt.title('Estimated Feedback/Feedforward Policies from R, D trials')
for sub in subjects:
    plt.loglog(freq[foIX],data[sub]['BDabs'],'b')
    plt.loglog(freq[foIX],data[sub]['FRabs'],'r')
plt.loglog(freq[foIX],np.absolute(np.divide(1,foM)),'r--',label='Model Inverse Prediction')
plt.ylabel('Gain')

plt.subplot(223)
for sub in subjects:
    plt.semilogx(freq[foIX],data[sub]['BDangle'],'b')
    plt.semilogx(freq[foIX],data[sub]['FRangle'],'r')
plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,'r--',label='Model Inverse Prediction')
plt.ylabel('Phase - deg')
plt.xlabel('Freq - Hz')


plt.subplot(222)
plt.title('Estimated Feedback/Feedforward Policies from R+D trials')
for sub in subjects:
    plt.loglog(freq[foIX],data[sub]['BRDabs'],'b')
    plt.loglog(freq[foIX],data[sub]['FRDabs'],'r')
plt.loglog(freq[foIX],np.absolute(np.divide(1,foM)),'r--',label='Model Inverse Prediction')
plt.ylabel('Gain')

plt.subplot(224)
for sub in subjects:
    plt.semilogx(freq[foIX],data[sub]['BRDangle'],'b')
    plt.semilogx(freq[foIX],data[sub]['FRDangle'],'r')
plt.semilogx(freq[foIX],np.angle(np.divide(1,foM))*180/m.pi,'r--',label='Model Inverse Prediction')
plt.ylabel('Phase - deg')
plt.xlabel('Freq - Hz')


# COHERENCE
datacoherence = {}
for sub in subjects:
    x =
    data[sub]['coherence'] = coherence(data[])
    data[sub]['BD'],data[sub]['FR'],data[sub]['BRD'],data[sub]['FRD'] = takeAverage(sub)
    data[sub]['BDabs'] = np.absolute(data[sub]['BD'])
    data[sub]['FRabs'] = np.absolute(data[sub]['FR'])
    data[sub]['BRDabs'] = np.absolute(data[sub]['BRD'])
    data[sub]['FRDabs'] = np.absolute(data[sub]['FRD'])

    data[sub]['BDangle'] = np.angle(data[sub]['BD'])*180/m.pi
    data[sub]['FRangle'] = np.angle(data[sub]['FR'])*180/m.pi
    data[sub]['BRDangle'] = np.angle(data[sub]['BRD'])*180/m.pi
    data[sub]['FRDangle'] = np.angle(data[sub]['FRD'])*180/m.pi

    dataAvg['BDavg'][:,numcount]=data[sub]['BDabs']
    dataAvg['FRavg'][:,numcount]=data[sub]['FRabs']
    dataAvg['BRDavg'][:,numcount]=data[sub]['BRDabs']
    dataAvg['FRDavg'][:,numcount]=data[sub]['FRDabs']

    angleAvg['BDavg'][:,numcount]=np.angle(data[sub]['BD'])*180/m.pi
    angleAvg['FRavg'][:,numcount]=np.angle(data[sub]['FR'])*180/m.pi
    angleAvg['BRDavg'][:,numcount]=np.angle(data[sub]['BRD'])*180/m.pi
    angleAvg['FRDavg'][:,numcount]=np.angle(data[sub]['FRD'])*180/m.pi
    numcount +=1

plt.show()
