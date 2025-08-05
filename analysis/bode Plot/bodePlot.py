# trying to make a bode plot
# http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=SystemAnalysis

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import math as m

# first order system
s1 = signal.lti([1],[1,1])
w,mag,phase = signal.bode(s1)
plt.figure()
plt.semilogx(w,mag)
plt.title('Standard feedback - magnitude')
plt.ylabel('Magnitude (dB)')
plt.xlabel('Frequency (Hz)')

plt.figure()
plt.semilogx(w,phase/360*2*m.pi)
plt.title('Standard feedback - phase')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')
plt.show()

# first order system
s1 = signal.lti([-1],[1,0])
w,mag,phase = signal.bode(s1)
plt.figure()
plt.semilogx(w,mag)
plt.title('M first order neg - magnitude')
plt.ylabel('Magnitude (dB)')
plt.xlabel('Frequency (Hz)')

plt.figure()
plt.semilogx(w,phase/360*2*m.pi)
plt.title('M first order neg - phase')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')

plt.figure()
s2 = signal.lti([-1,0],[1])
w,mag,phase = signal.bode(s2)
plt.semilogx(w,mag)
plt.title('M inv neg - magnitude')
plt.ylabel('Magnitude (dB)')
plt.xlabel('Frequency (Hz)')

plt.figure()
plt.semilogx(w,phase/360*2*m.pi)
plt.title('M inv neg - phase')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')

# first order system
s1 = signal.lti([1],[1,0])
w,mag,phase = signal.bode(s1)
plt.figure()
plt.semilogx(w,mag)
plt.title('M first order (pos) - magnitude')
plt.ylabel('Magnitude (dB)')
plt.xlabel('Frequency (Hz)')

plt.figure()
plt.semilogx(w,phase/360*2*m.pi)
plt.title('M first order pos - phase')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')

plt.figure()
s2 = signal.lti([1,0],[1])
w,mag,phase = signal.bode(s2)
plt.semilogx(w,mag)
plt.title('M inv pos - magnitude')
plt.ylabel('Magnitude (dB)')
plt.xlabel('Frequency (Hz)')

plt.figure()
plt.semilogx(w,phase/360*2*m.pi)
plt.title('M inv - phase')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')
plt.show()
#
# # second order system
# s2 = signal.lti([1],[1,1,1])
# w,mag,phase = signal.bode(s2)
# plt.figure()
# plt.semilogx(w,mag)
# plt.ylabel('magnitude (dB)')
# plt.title('Second order system - magnitude')
#
# #plt.figure()
# plt.semilogx(w,phase)
# plt.title('Second order system - phase')
# plt.xlabel('Frequency (rad)')
# plt.ylabel('Phase (rad)')
# plt.show()
#
# # plot first order system model simulation
# Kp = 1.0
# taup = 1.0
#
# # transfer function
# num = [Kp]
# den = [taup,1]
# sys1 = signal.TransferFunction(num,den)
# t1,y1 = signal.step(sys1)
#
# # state space
# A = -1.0/taup
# B = Kp/taup
# C = 1.0
# D = 0.0
# sys2 = signal.StateSpace(A,B,C,D)
# t2,y2 = signal.step(sys2)
#
#  # ODE integrator
# def model3(y,t):
#      u = 1
#      return(-y + Kp * u)/taup
# t3 = np.linspace(0,14,100)
# y3 = odeint(model3,0,t3)
#
# plt.figure(1)
# plt.plot(t1,y1,'b--',linewidth=3,label='Transfer Fcn')
# plt.plot(t2,y2,'g:',linewidth=2,label='State Space')
# plt.plot(t3,y3,'r-',linewidth=1,label='ODE Integrator')
# plt.xlabel('Time')
# plt.ylabel('Response (y)')
# plt.legend(loc='best')
plt.show()
