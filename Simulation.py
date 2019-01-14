import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp

alpha = 0.4*np.pi
beta = 0.01*np.pi
eps = 0.001
N = 3
state = [np.array([0,1,2]),np.array([0,0,0,0,0,0,0,0,0])]

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

def derivs(state,t):
	Phis,Kappas = state
	DelPhis = (Phis-Phis[:,np.newaxis]).T
	DKappa = -eps*(np.sin(DelPhis+beta)+Kappas)
	DelPhis += alpha
	DelPhis = Kappas*np.sin(DelPhis)
	DPhis = DelPhis.sum(axis=1)/Phis.shape[0]
	derivatives = [DPhis,DKappa]
	return derivatives

def Calccoords(state):
	while state[0]>2*np.pi:
		state[0]-=2*np.pi
	while state[1]>2*np.pi:
		state[1]-=2*np.pi
	while state[2]>2*np.pi:
		state[2]-=2*np.pi
	state.sort()
	s0 = 0
	s1 = state[1]-state[0]
	s2 = state[2]-state[0]
	l3 = (2*np.pi-s2)/(2*np.pi)
	l2 = np.abs(s2-s1)/(2*np.pi)
	l1 = np.abs(s1-s0)/(2*np.pi)
	Sum = l1+l2+l3
	return [l1,l2,l3]

def Coordstopointinspace(l):
	x = (l[1]+l[0]*0.5)/(l[0]+l[1]+l[2])
	y = l[0]*np.sqrt(3)/2/(l[0]+l[1]+l[2])
	return [x,y]

if __name__=="__main__":
	Phi0 = np.arange(3,dtype='float64')
	Phi0 *= 2./3*np.pi
	Kappa0 = np.ones((3,3),dtype='float64')
	Kappa0 *= np.sin(beta)
	state0 = [[Phi0],[Kappa0]]
#	state0 = np.zeros((12,))
#	state0 -= np.array([0, -2./3*np.pi, -4./3*np.pi, np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta), np.sin(beta)])
#	t = np.arange(0,1000,0.1)
	print (Gamma(1./3.)-alpha-beta)/np.pi
#	states = odeint(derivs,state0,t)
	states = state0
	t=0
	while t<1e6:
		calc_state = [states[0][-1],states[1][-1]]
#		print calc_state
		
		d = derivs(calc_state,0)
#		print d
		states[0].append(np.array(states[0][-1])+0.001*d[0])
		states[1].append(np.array(states[1][-1])+0.001*d[1])
		t+=1
		if t%1e3:
			print t/1000.
	Phis = [[],[],[]]
	for i in xrange(len(Phis)):
		for j in xrange(len(states[0])):
			Phis[i].append(states[0][j][i])
			
	for Set in Phis:
		plt.plot(Set)
	plt.show()
	"""
	Xvals = []
	Yvals = []
	for i in xrange(len(Phis[0])/1):
		state = [Phis[0][i],Phis[1][i],Phis[2][i]]
		Coords = Calccoords(state)
		point = Coordstopointinspace(Coords)
		#print Coords
		Xvals.append(point[0])
		Yvals.append(point[1])
	
	plt.scatter(Xvals,Yvals, s=0.05)
	plt.plot([0,1],[0,0],lw=0.2,c="k")
	plt.plot([0,0.5],[0,np.sqrt(3)/2],lw=0.2,c="k")
	plt.plot([0.5,1],[np.sqrt(3)/2,0],lw = 0.2, c="k")
	plt.show()
	"""
