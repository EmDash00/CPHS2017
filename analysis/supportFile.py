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
from globalVars import primes,IX,primesE,primesO,IXE,IXO,freq,tempfreq,TEMPFREQ,s,M
import sys

# go through list of filenames and find FFT of each timestream of data
def findFFT(listofIDs,IX,subject):
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
    import sys
    from scipy.optimize import curve_fit
    from findFFT import findFFT
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}

    matplotlib.rc('font', **font)

    # Define
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
    M = 1.0/(np.multiply(s,s)-s) # plus or minus s???
    M = np.asarray(M)
    #sys.path.insert(0,'C:\\Users\\BRL\\Momona\'s Google Drive\\Yamagami Lab_\\Slider Project\\hcps\\develop\\')
    trials = dict()
    for id in listofIDs:

      fis = sorted(glob.glob(os.path.join('D:\\Momona\\Google Drive\\NEW Yamagami Lab\\hcps\\experiment\\data\\',subject,id+'.npz')))#'*_'+protocol+'_'+id+'.npz')))
      if len(fis) > 1:
        dbg('WARNING -- repeated trials for id ='+id)
      assert len(fis) > 0, 'ERROR -- no data for id ='+id
      fi = fis[-1]
      #print('LOAD '+fi)
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

# return time and frequency domain measures
takeAvgBA = True
def plotThingsEO(M,Bothid,Evenid,Oddid,IX,IXE,IXO,subject):
    #global err
    err = np.zeros((int(len(Bothid)-1),),dtype=complex)
    #global timeerr
    timeerr = np.zeros((int(len(Bothid)),))

    # find ffts
    both,BOTH = findFFT(Bothid,IX,subject)

    for num,id in enumerate(Bothid):
        timeerr[num] = np.mean((np.asarray(both[id]['outs'])-np.asarray(both[id]['refs']))**2)
    BA = np.zeros((int(len(Bothid)-1),len(IX)),dtype=np.complex_) # 29 trials by 8 freqs
    FA = np.zeros((int(len(Bothid)-1),len(IX)),dtype=np.complex_)
    FAavg = np.zeros((int(len(Bothid)),len(IX)),dtype=np.complex_)

    for num,id in enumerate(Bothid[:-1]):
        firstID = id
        secID = Bothid[num+1]
        Hdu = np.zeros(BOTH[firstID]['Hdu'].shape,dtype=np.complex_)
        Hru = np.zeros(Hdu.shape,dtype=np.complex_)
        if firstID in Evenid:
#         if sum(abs(BOTH[firstID]['REFS'][IXE])) < 1e-12: # this means that the first ID is an ODD trial
            Hru[::2] = BOTH[firstID]['Hru'][::2]
            Hdu[1::2] = BOTH[firstID]['Hdu'][1::2]
            Hru[1::2] = BOTH[secID]['Hru'][1::2]
            Hdu[::2] = BOTH[secID]['Hdu'][::2]
        elif firstID in Oddid:#sum(abs(BOTH[secID]['REFS'][IXE])) < 1e-12: # this means that the first ID is an EVEN trial
            Hru[::2] = BOTH[secID]['Hru'][::2]
            Hdu[1::2] = BOTH[secID]['Hdu'][1::2]
            Hru[1::2] = BOTH[firstID]['Hru'][1::2]
            Hdu[::2] = BOTH[firstID]['Hdu'][::2]
        else:
            print('error: something is wrong')
        BA[num,:] = np.divide(-Hdu,np.multiply(M,np.ones(Hdu.shape,dtype=np.complex_)+Hdu))
        FA[num,:] = np.multiply(Hru,(1+0j)+np.multiply(BA[num,:],M))-BA[num,:]

    for idNum in range(0,len(err),1):
        err[idNum] = np.mean(((abs(abs(FA[idNum])-abs(np.divide(np.ones(M.shape),M))))**2))

    return err,timeerr,BA,FA

# find files for specific subject
def findFilename(subject):
    fis = glob.glob('D:\\Momona\\Google Drive\\NEW Yamagami Lab\\hcps\\experiment\\data\\'+subject+'\\*.npz')
    #fis = glob.glob('C:\\Users\\BRL\\Momona\'s Google Drive\\Yamagami Lab_\\Slider Project\\hcps\\develop\\data/'+subject+'\\*.npz')
    #print(fis)
    fis = [os.path.basename(fi) for fi in fis]

    ids = sorted(list(set([fi.strip('.npz') for fi in fis])))
    #ids = [id for id in ids if ('.csv') not in id[-4:]]
    ids = [id for id in ids if ('rst2') not in id]
    ids = [id for id in ids if ('rst1') not in id]
    ids = [id for id in ids if ('rst0') not in id]
    ids = [id for id in ids if ('react') not in id]
    #ids = [id for id in ids if ('st2', 'st1','rej','oob','rst','max','st1','st2') not in id]
    #ids = [id for id in ids if ('so') in id]
    Aid = [id for id in ids if ('r-sos+A_d-sos+A') in id]
    Evenid = [id for id in ids if ('r-sos+E') in id or ('r-sos-E') in id]
    Oddid = [id for id in ids if ('r-sos-O') in id or ('r-sos+O') in id]
    Bothid = [id for id in ids if ('sos+A') not in id]
    return ids,Evenid,Oddid

# combine all functions to get variables of interest
def geterr(subject,M):
    ids,Evenid,Oddid = findFilename(subject)
    err, timeerr,BA,FA = plotThingsEO(M,ids,Evenid,Oddid,IX,IXE,IXO,subject)
    return err, timeerr,BA,FA
