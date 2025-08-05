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

      fis = sorted(glob.glob(os.path.join('C:\\Users\\BRL\\Momona\'s Google Drive\\Yamagami Lab_\\Slider Project\\hcps\\experiment\\data\\',subject,id+'.npz')))#'*_'+protocol+'_'+id+'.npz')))
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
