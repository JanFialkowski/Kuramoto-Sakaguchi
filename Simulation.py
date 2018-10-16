import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

alpha = 0
beta = -np.pi
eps = 0.001
N = 3
state = [0,1,2,0,0,0,0,0,0,0,0,0]

def deltkappa(state,i,j):
	return -eps*(np.sin(state[i-1]-state[j-1]+beta)+state[3*i+j-1])

def deltphi(state,i):
	Res = 0
	for j in xrange(N):
		Res += state[3*i+j-1]*np.sin(state[i-1]-state[j-1]+alpha)
	Res /= N
	return Res

def f(state,t):
	derivs = []
	for i in xrange(N):
		derivs.append(deltphi(state,i))
	for i in xrange(N):
		for j in xrange(N):
			derivs.append(deltkappa(state,i,j))
	return derivs

state0 = np.zeros((12,))
state0 -= np.array([0,0,0,np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta)])
t = np.arange(0,1000,0.001)

states = odeint(f,state0,t)
Phis = [[],[],[]]
for i in xrange(len(Phis)):
	for j in xrange(len(states)):
		Phis[i].append(states[j][i]/np.pi)
		
for Set in Phis:
	plt.plot(t,Set)
plt.show()
