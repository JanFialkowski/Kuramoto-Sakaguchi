import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

alpha = 0.01*np.pi
beta = 0.3*np.pi
eps = 0.1
N = 3
state = [0,1,2,0,0,0,0,0,0,0,0,0]

def Gamma(n):
	delta = 2*(alpha+beta)
	s = np.sin(delta)
	c = np.cos(delta)
	G = n*s
	G2 =((1+c)*n-1)
	Out=np.arctan2(G,G2)
	return Out

def deltkappa(state,i,j): #i,j are oscillatorindices in [1,N]
	return -eps*(np.sin(state[i-1]-state[j-1]+beta)+state[3*i+j-1])

def deltphi(state,i):#i in [1,N]
	Res = 0
	for j in xrange(N):#j in [0,N-1]
		Res += state[3*i+j]*np.sin(state[i-1]-state[j]+alpha)
	Res /= -N
	return Res

def f(state,t):
	derivs = []
	for i in xrange(N):
		derivs.append(deltphi(state,i+1))#i in [0,N-1]
	for i in xrange(N):
		for j in xrange(N):
			derivs.append(deltkappa(state,i+1,j+1))#i,j in [0,N-1]
	return derivs

def Calccoords(state):
	l1 = (state[2]-state[1])/np.pi%1
	l2 = (state[3]-state[2])/np.pi%1
	return [l1,l2]


state0 = np.zeros((12,))
state0 -= np.array([0,0,0.1,np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta),np.sin(beta)])
t = np.arange(0,10000,0.01)

print (Gamma(1./3.)-alpha-beta)/np.pi
states = odeint(f,state0,t)
Phis = [[],[],[]]
for i in xrange(len(Phis)):
	for j in xrange(len(states)):
		Phis[i].append(states[j][i]/np.pi)
		
for Set in Phis:
	plt.plot(t,Set)
plt.show()
