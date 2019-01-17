import numpy as np

def Gamma(n, alpha, beta):
	delta = 2*(alpha+beta)
	s = np.sin(delta)
	c = np.cos(delta)
	G = n*s
	G2 =-((1+c)*n-1)
	Out=np.arctan2(G,G2)
	return Out

def derivs(state,t, alpha, beta, eps):
	N = int(-0.5+np.sqrt(0.25+len(state)))
	Phis = state[:N]
	Kappas = state[N:].reshape(N,N)
	DelPhis = (Phis-Phis[:,np.newaxis]).T
	DKappa = -eps*(np.sin(DelPhis+beta)+Kappas)
	DelPhis += alpha
	DelPhis = Kappas*np.sin(DelPhis)
	DPhis = DelPhis.sum(axis=1)/Phis.shape[0]*(-1)
	derivatives = np.concatenate((DPhis,DKappa.flatten()))
	if t%100==0:
		print t/100.
	return derivatives

def L1(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(Theta)-n*np.sin(-2*Psi+Theta)+2*eps)
	q = -eps*((1-n)*np.sin(Theta)+n*np.sin(-2*Psi+Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)

def L2(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(2*Psi+Theta)-n*np.sin(Theta)+2*eps)
	q = -eps*((1-n)*np.sin(2*Psi+Theta)+n*np.sin(Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)

def L3(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(2*Psi+Theta)-n*np.sin(-2*Psi+Theta)+2*eps)
	q = -eps*((1-n)*np.sin(2*Psi+Theta)+n*np.sin(-2*Psi+Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)

if __name__ == "__main__":
	alpha = np.linspace(0,np.pi/2)
	beta = np.linspace(-np.pi,np.pi)
	eps = 0.01
	n= 0.1
	print L1(n,alpha,beta,eps)
