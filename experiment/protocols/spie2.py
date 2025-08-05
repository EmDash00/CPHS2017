import numpy as np
import pylab as plt
import matplotlib 
from matplotlib import rc

from dynamics import fo, so, zd11, zd12

from references import spline as ref_
from references import spline_interp

from disturbances import zero as dis

nls = dict()
# sigmoidal nonlinearity (maps [-.5,.5] to [-.5,.5])
nls['sig'] = lambda h,s=4. : np.arctan(s*h)/np.arctan(s/2.)/2.
# inverse of sigmoidal nonlinearity -- ln \circ nl = \id
nls['gis'] = lambda y,s=4. : np.tan(2*np.arctan(s/2.)*y)/s
# vector fields
vfs = ['fo']
# reference goes from -jump to jump for each jump in jump
jumps = [0.1,0.2]
# number of replicates
I = 10 
# time to complete task
dt = 1.
#
scale = dict()
scale['fo'] = spline_interp(.5,0.,1.,-max(jumps),+max(jumps),0.,0.,return_dy=True)[1]
scale['so'] = 0.25*spline_interp(0,0.,1.,-max(jumps),+max(jumps),0.,0.,return_ddy=True)[1]
scale['zd11'] = 0.25*spline_interp(0,0.,1.,-max(jumps),+max(jumps),0.,0.,return_ddy=True)[1]
scale['zd12'] = 0.25*spline_interp(0,0.,1.,-max(jumps),+max(jumps),0.,0.,return_ddy=True)[1]
#
def trial_gen(subject,protofile):
  for vf in vfs:
    for nl in ['gis','sig']:
      ln = nl[::-1] # inverse nonlinearity
      # apply nonlinearity to reference
      ref = lambda t,trial : nls[nl](ref_(t,trial))
      for jump in jumps:
        for i in range(I):
          for s in [+1,-1]:
            id = '%s_%s_j%+0.1f_i%d' % (vf,nl,s*jump,i)
            if vf == 'fo':
              state = [-s*jump]
              out = lambda x : nls[nl](x[0])
            elif vf == 'so':
              state = [-s*jump,0.]
              out = lambda x : nls[nl](x[0])
            elif vf == 'zd11':
              state = [-s*jump,-s*jump]
              out = lambda x : nls[nl](.5*(x[0]+x[1]))
            elif vf == 'zd12':
              state = [-s*jump,-s*jump,0.]
              out = lambda x : nls[nl](.5*(x[0]+x[1]))

            #
            pts = np.vstack(([-np.inf,1.*dt, -s*jump,-s*jump],
                             [1.*dt,  2.*dt, -s*jump,+s*jump],
                             [2.*dt,  np.inf,+s*jump,+s*jump],
                           ))
            #
            trial = dict(id=id,pts=pts,init=state,vf=vf,ref=ref,dis=dis,out=out,
                         scale=scale[vf])
            yield trial

# plotting globals
lw = 1
col_out = 'indigo'
col_ref = 'darkorange'
col_p_gain = 'dodgerblue'
col_pd_gain = 'dodgerblue'
col_inp = col_out
alpha_trials = 0.0
alpha_fill = 0.2
#
xlim = (1.,3.)
#ylim_out = (-.3,+.3)
ylim_out = (-0.1,+.5)
#
prc = 25 
prcs = [50,prc/2,100-prc/2] 
do_bootstrap_means = False
#
# multiplicative factor for "aggressiveness" of pure feedback
aggro = 1.
do_gains = True
#do_gains = False
#
# pure feedback simulation
def gain_sim(gains,vf,nl,ln,time,state0,ref):
  state = np.nan*np.zeros((time.size,state0.size))
  state[0] = state0
  state[1] = state0
  inp = np.nan*time
  inp[0] = 0.
  inp[-1] = 0.
  for ti in range(1,time.size-1):
    dt = time[ti+1] - time[ti]
    t = time[ti]
    ede = np.vstack(( ref[ti] - nl(state[ti]),
                    ( (ref[ti] - ref[ti-1]) - (nl(state[ti]) - nl(state[ti-1])) ) / dt ))
    #ede = np.vstack(( ln(ref[ti]) - state[ti],
    #                ( (ln(ref[ti]) - ln(ref[ti-1])) - (state[ti] - state[ti-1]) ) / dt ))
    inp[ti] = np.dot( gains, ede )
    state[ti+1] = state[ti] + dt * vf(t,state[ti],inp[ti])
  #1/0
  return state,inp


def analysis(trials,title='',rescale=False,**util):
  #
  ids = sorted(trials.keys())
  #
  time_ = trials[ids[0]]['time_']
  samps = np.arange((time_ <= xlim[0]).nonzero()[0][-1], (time_ < xlim[1]).nonzero()[0][-1])
  times = time_[samps]
  # 
  figs = dict()
  axs = dict()
  lines = dict()
  #
  errors = dict(); derrors = dict(); inps = dict()
  errors['all'] = []; derrors['all'] = []; inps['all'] = []
  p_gains = dict(); p_gains['all'] = []
  pd_gains = dict(); pd_gains['all'] = []
  for vi,vf in enumerate(vfs):
    #
    # gather subjects' inputs and tracking errors from spie1
    for s in ['+','-','+/-']:
      for i,id in enumerate(ids):
        subj,proto,_ = id.split('_',2)
        if not proto == 'spie1':
          continue
        subj,proto,v,j,i = id.split('_')
        if not v == vf: # or not j[1] in s:
          continue
        trial = trials[id]
        time = trial['time_']
        # Eatai:  I think we need to restrict to end times
        #gain_samps = (time >= 3.0).nonzero()[0]
        gain_samps = ((trial['time_'] >= 1.5)*(trial['time_'] <= 2.5)).nonzero()[0]
        # collect output errors and inputs for subject
        error = trial['ref_'] - trial['out_']
        derror = np.hstack((0.,np.diff(error)/np.diff(time))) # + np.random.randn(time.size)*1e-2
        inp = trial['inp_']
        if subj not in errors:
          errors[subj] = []; derrors[subj] = []; inps[subj] = []
        #
        errors[subj].append(error[gain_samps])
        derrors[subj].append(derror[gain_samps])
        inps[subj].append(inp[gain_samps])
        #
        errors['all'].append(error[gain_samps])
        derrors['all'].append(derror[gain_samps])
        inps['all'].append(inp[gain_samps])
        # subtract feedforward contribution
        #if subj+'_' not in errors:
        #  errors[subj+'_'] = []; inps[subj+'_'] = []
        #inp[1:] -= np.diff(trial['ref_']) / np.diff(trial['time_'])
        #errors[subj+'_'].append(error[gain_samps])
        #inps[subj+'_'].append(inp[gain_samps])
    # estimate proportional feedback gain from errors and inputs
    for subj in errors.keys():
      # solve for gain via linear least-squares: 
      #   error = gain * inp
      error = np.hstack(errors[subj])
      derror = np.hstack(derrors[subj])
      inp = np.hstack(inps[subj])
      U = inp; E = error # np.vstack((error,derror))
      p_gains[subj] = np.dot( np.dot(U, E.T), 
                              1./np.dot(E,E.T) )
      U = inp; E = np.vstack((error,derror))
      pd_gains[subj] = np.dot( np.dot(U, E.T), 
                               np.linalg.inv(np.dot(E,E.T)) )
    open('p_gains.txt','w').write(str(p_gains))
    open('pd_gains.txt','w').write(str(pd_gains))
    print p_gains
    print pd_gains
    #print [(k,v) for k,v in p_gains.items() if '_' not in k]
    #print [(k,v) for k,v in pd_gains.items() if '_' not in k]
    #1/0
    #
    for ni,_nl in enumerate(['gis','sig']):
      inps,outs = dict(),dict()
      #
      #ylim_inp = scale[vf] * (3./2.) * np.array([-0.75,+0.75])
      #ylim_inp = scale[vf] * (3./2.) * np.array([-1.,+1.])
      ylim_inp = scale[vf] * (3./2.) * np.array([-.25,+1.])
      #
      figs[vf+'_'+_nl] = plt.figure(2+vi*len(vfs)+ni,figsize=(12,6)); plt.clf()
      axs[vf] = dict()
      axs[vf]['+'] = [plt.subplot(2,3,1),plt.subplot(2,3,4)]
      axs[vf]['-'] = [plt.subplot(2,3,2),plt.subplot(2,3,5)]
      axs[vf]['+/-'] = [plt.subplot(2,3,3),plt.subplot(2,3,6)]
      #
      for s in axs[vf].keys():
        for i,id in enumerate(ids):
          subj,proto,_ = id.split('_',2)
          if not proto == 'spie2':
            continue
          subj,proto,v,nl,j,i = id.split('_')
          ln = nl[::-1] # inverse nonlinearity
          if not v == vf or not j[1] in s or not nl == _nl:
            continue
          sgn = float(j[1]+'1')
          jump = float(j[2:])
          sc = 1.
          if rescale:
            sc = max(jumps)/jump
          trial = trials[id]
          ax0,ax1 = axs[vf][s]
          #
          ylabel_inp = ''; ylabel_out = ''
          if s == '+':
            ylabel_inp = r'input ($u$)'
            ylabel_out = r'output ($y$)'
          # compute pure feedback input
          if vf == 'fo' and do_gains:
            p_gain_state,p_gain_inp = gain_sim(
                                           [aggro*p_gains['all'],0.],
                                           eval(vf),
                                           nls[nl],
                                           nls[ln],
                                           trial['time_'],
                                           trial['state_'][0],
                                           trial['ref_'])
            pd_gain_state,pd_gain_inp = gain_sim(
                                           aggro*pd_gains['all'],
                                           eval(vf),
                                           nls[nl],
                                           nls[ln],
                                           trial['time_'],
                                           trial['state_'][0],
                                           trial['ref_'])
          #1/0
          # out
          #ax0.plot(trial['time_'],sc*2*jump,'k--',lw=lw,zorder=-10,alpha=alpha_trials)
          if do_gains:
            lines['p_out']  = ax0.plot(trial['time_'],sgn*sc*(p_gain_state)+jump,'-',lw=2*lw,color=col_p_gain,zorder=-10)
            lines['pd_out'] = ax0.plot(trial['time_'],sgn*sc*(pd_gain_state)+jump,'--',lw=2*lw,color=col_pd_gain,zorder=-10)
          lines['ref'] = ax0.plot(trial['time_'],sgn*sc*nls[ln](trial['ref_'])+jump,lw=2*lw,color=col_ref,zorder=-10)
          lines['out'] = ax0.plot(trial['time_'],sgn*sc*nls[ln](trial['out_'])+jump,'-',lw=lw,color=col_out,zorder=0,alpha=alpha_trials)
          if s == '+':
            xlabel_out = 'up'
            yticklabels = None
          elif s == '-':
            xlabel_out = 'down'
            yticklabels = []
          else:
            xlabel_out = 'pooled up/down'
            yticklabels = []
          util['format_axes'](ax0,title=u''+title+' '+vf+' '+s,
                              xlim=xlim,ylim=ylim_out,
                              xlabel=xlabel_out,ylabel=ylabel_out,
                              xticklabels=[],yticklabels=yticklabels)
          if j[1:] not in outs:
            outs[j[1:]] = []
          outs[j[1:]].append(sgn*sc*nls[ln](trial['out_'][samps])+jump)
          # inp
          if vf == 'fo' or vf == 'zd12':
            # pure feedforward input
            pred = np.diff(sgn*sc*nls[ln](trial['ref_'])) / np.diff(trial['time_'])
            pred = np.hstack((pred,[0.]))
            #
          if vf == 'so':
            pred = np.diff(np.diff(sgn*sc*nls[ln](trial['ref_']))) / np.diff(trial['time_'])[1:]**2
            pred = np.hstack((pred,[0.,0.]))
          if do_gains:
            lines['p_inp']  = ax1.plot(trial['time_'],sgn*sc*p_gain_inp,'-',lw=2*lw,color=col_p_gain,zorder=-10)
            lines['pd_inp'] = ax1.plot(trial['time_'],sgn*sc*pd_gain_inp,'--',lw=2*lw,color=col_pd_gain,zorder=-10)
          lines['pred'] = ax1.plot(trial['time_'],pred,'-',lw=2*lw,color=col_ref,zorder=-10)
          lines['inp']  = ax1.plot(trial['time_'],sgn*sc*trial['inp_'],'-',lw=lw,color=col_out,zorder=0,alpha=alpha_trials)
          if s == '+':
            xlabel_inp = ''
            yticklabels = None
          elif s == '-':
            xlabel_inp = 'time (sec)'
            yticklabels = []
          else:
            xlabel_inp = ''
            yticklabels = []
          util['format_axes'](ax1,xlabel=xlabel_inp,
                              xlim=xlim,ylim=ylim_inp,
                              ylabel=ylabel_inp,
                              yticklabels=yticklabels)
          if j[1:] not in inps:
            inps[j[1:]] = []
          inps[j[1:]].append(sgn*sc*trial['inp_'][samps])
          #
        # distributions
        ax0,ax1 = axs[vf][s]
        for jump in jumps:
          outs_ = []; inps_ = []
          for sgn in ['+','-']:
            if sgn in s:
              outs_.append(outs[sgn+str(jump)])
              inps_.append(inps[sgn+str(jump)])
          outs_ = np.vstack(outs_)
          inps_ = np.vstack(inps_)
          if len(outs_) == 0:
            continue
          #
          if do_bootstrap_means:
            bootstrap_means = util['bootstrap_mean'](outs_)
            out = np.percentile(bootstrap_means,prcs,axis=0)
            #bootstrap_trials = util['bootstrap_trials'](outs[str(jump)])
            #out = np.percentile(bootstrap_trials,[50,1,99],axis=0)
          else:
            out = np.percentile(outs_,prcs,axis=0)
          lines['out_fill'] = ax0.fill_between(times,out[1],out[2],color=col_out,zorder=10,alpha=alpha_fill)
          lines['out_bd']   = ax0.plot(times,out[1:].T,'--',lw=1*lw,color=col_out,zorder=10)
          lines['out_mean'] = ax0.plot(times,out[0],   '-', lw=2*lw,color=col_out,zorder=10)
          #
          if do_bootstrap_means:
            bootstrap_means = util['bootstrap_mean'](inps_)
            inp = np.percentile(bootstrap_means,prcs,axis=0)
          else:
            inp = np.percentile(inps_,prcs,axis=0)
          lines['inp_fill'] = ax1.fill_between(times,inp[1],inp[2],color=col_inp,zorder=10,alpha=alpha_fill)
          lines['inp_bd']   = ax1.plot(times,inp[1:].T,'--',lw=1*lw,color=col_inp,zorder=10)
          lines['inp_mean'] = ax1.plot(times,inp[0],   '-', lw=2*lw,color=col_inp,zorder=10)
        # legend
        if s == '-':
          lines['ref'][0].set_label('reference')
          lines['out_mean'][0].set_label('operator average')
          #lines['out_fill'].set_label('operator 50%')
          lines['out_bd'][0].set_label('operator quartiles')
          lines['p_inp'][0].set_label('P feedback')
          lines['pd_inp'][0].set_label('PD feedback')
          ls = (lines['ref'][0],lines['out_mean'][0],lines['out_bd'][0],lines['p_inp'][0],lines['pd_inp'][0])
          lb = ['reference','operator average','operator quartiles','P feedback','PD feedback']
          leg_pos = (0.5,1.05)
          ax0.legend(ls,lb,ncol=6,loc='lower center',bbox_to_anchor=leg_pos)
  
  figs['nl'] = plt.figure(2+vi*len(vfs)+ni+1,figsize=(6,6))
  plt.clf()
  ax = plt.subplot(1,1,1)
  ax.grid('on'); ax.axis('equal')
  lim = (-.5,.5)
  y = np.linspace(*lim)
  ax.plot(y,y,'k--',lw=2*lw,label=u'$h(y) = y$')
  ax.plot(y,nls['sig'](y),'b',lw=2*lw,label=u'$h(y) = \operatorname{sig}(y)$')
  ax.plot(y,nls['gis'](y),'r',lw=2*lw,label=u'$h(y) = \operatorname{gis}(y)$')
  ax.legend(loc='upper left')
  util['format_axes'](ax,
                      xlim=lim,ylim=lim,
                      xlabel=u'$y$',ylabel=u'$h(y)$')
 
  #
  return figs

