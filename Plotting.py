import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as manimation
import Kuramoto as K
from mpl_toolkits.mplot3d import Axes3D
import Simulation as S
#from scipy.integrate import solve_ivp

alpha = 0.4*np.pi
beta = -0.15*np.pi
eps = 0.01
N=10
N2 = 1.

Phi0 = np.zeros(N)
G = S.Gamma(N2/N)+alpha+beta
#Phi0[0]= -G
#Phi0[1] = np.pi
Kappa0 = -np.sin((Phi0-Phi0[:,np.newaxis]).T+beta).flatten()
Phi0[0] += 0.001
state0 = np.concatenate((Phi0,Kappa0))
t = np.arange(0,3000,.1)
print (-G)/np.pi
Phis, Kappa = S.Integrator(state0, t, 10)
state0[0]-=0.002
Phis2, Kappa2 = S.Integrator(state0, t, 10)
state0[0]+=0.001
#"""
Phi0 = np.zeros(N)#np.random.rand(50)*2*np.pi*0.01
G = S.Gamma(N2/N)+alpha+beta
Phi0[0]= np.pi
#Phi0[1] = np.pi
Kappa0 = -np.sin((Phi0-Phi0[:,np.newaxis]).T+beta).flatten()
Phi0[0] += 0.001
state0 = np.concatenate((Phi0,Kappa0))
print (-G)/np.pi
Phis3, Kappa3 = S.Integrator(state0, t, 10)
state0[0]-=0.002
Phis4, Kappa4 = S.Integrator(state0, t, 10)
state0[0]+=0.001
#"""
fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")
#ax2 = fig.add_subplot(122,projection = "3d")
Theta = (np.array(Phis[0])-np.array(Phis[1]))/np.pi
Theta2 = (np.array(Phis2[0])-np.array(Phis2[1]))/np.pi
Theta3 = (np.array(Phis3[0])-np.array(Phis3[1]))/np.pi
Theta4 = (np.array(Phis4[0])-np.array(Phis4[1]))/np.pi
Ind1 = np.where(Theta3<=2)[0].tolist()
Ind2 = np.where(Theta3>2)[0].tolist()
Ind3 = np.where(Theta2<=0)[0].tolist()
Ind4 = np.where(Theta2>0)[0].tolist()
print len(Kappa3)
ax.plot(Kappa[0], Kappa[1], Theta, c="blue")
ax.plot(xs=np.array(Kappa2)[0,:Ind3[0]], ys=np.array(Kappa2)[1,:Ind3[0]], zs=Theta2[:Ind3[0]], c="brown")
ax.plot(xs=np.array(Kappa2)[0,Ind3[0]:], ys=np.array(Kappa2)[1,Ind3[0]:], zs=Theta2[Ind3[0]:], c="brown")
ax.plot(xs=np.array(Kappa3)[0,Ind1], ys=np.array(Kappa3)[1,Ind1], zs=Theta3[Ind1], c="green")
ax.plot(xs=np.array(Kappa3)[0,Ind2], ys=np.array(Kappa3)[1,Ind2], zs=Theta3[Ind2]-2, c="green")
ax.plot(Kappa4[0], Kappa4[1], Theta4, c="red")
ax.scatter(-np.sin(-G-Phi0[1]+np.pi+beta), -np.sin(Phi0[1]+G-np.pi+beta), (-G-Phi0[1])/np.pi+1, c="black")
ax.scatter(-np.sin(-G-Phi0[1]+beta), -np.sin(Phi0[1]+G+beta), (-G-Phi0[1])/np.pi, c="black")
ax.scatter(-np.sin(-G-Phi0[1]+beta), -np.sin(Phi0[1]+G+beta), (-G-Phi0[1])/np.pi+2, c="black")
#ax.scatter(-np.sin(beta),-np.sin(beta),0,c="orange")
ax.scatter(-np.sin(beta+np.pi),-np.sin(beta+np.pi),1,c="orange")
ax.scatter(-np.sin(beta),-np.sin(beta),0,c="orange")
ax.set_xlabel("Kappa12")
ax.set_ylabel("Kappa21")
ax.set_zlabel("Theta")
#ax2.set_xlabel("Kappa12")
#ax2.set_ylabel("Kappa21")
#ax2.set_zlabel("Theta")
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.show()
fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")
ax.plot(xs=np.array(Kappa3)[0,Ind1], ys=np.array(Kappa3)[1,Ind1], zs=Theta3[Ind1], c="green")
ax.plot(xs=np.array(Kappa3)[0,Ind2], ys=np.array(Kappa3)[1,Ind2], zs=Theta3[Ind2]-2, c="green")
ax.plot(Kappa4[0], Kappa4[1], Theta4, c="red")
ax.scatter(-np.sin(-G-Phi0[1]+np.pi+beta), -np.sin(Phi0[1]+G-np.pi+beta), (-G-Phi0[1])/np.pi+1, c="black")
ax.scatter(-np.sin(-G-Phi0[1]+beta), -np.sin(Phi0[1]+G+beta), (-G-Phi0[1])/np.pi+2, c="black")
#ax.scatter(-np.sin(beta),-np.sin(beta),0,c="orange")
ax.scatter(-np.sin(beta+np.pi),-np.sin(beta+np.pi),1,c="orange")
ax.scatter(-np.sin(beta),-np.sin(beta),0,c="orange")
ax.set_xlabel("Kappa12")
ax.set_ylabel("Kappa21")
ax.set_zlabel("Theta")
plt.show()
R1 = np.abs(S.CalcR2(np.array(Phis)))
R2 = np.abs(S.CalcR2(np.array(Phis2)))
R3 = np.abs(S.CalcR2(np.array(Phis3)))
R4 = np.abs(S.CalcR2(np.array(Phis4)))
Zeta = [np.array(Phis[0])/np.pi-np.array(Phis[1])/np.pi,np.array(Phis[0])/np.pi-np.array(Phis[2])/np.pi, "blue"]
Zeta4 = [np.array(Phis4[0])/np.pi-np.array(Phis4[1])/np.pi,np.array(Phis4[0])/np.pi-np.array(Phis4[2])/np.pi,"red"]
Zeta2 = [np.array(Phis2[0])/np.pi-np.array(Phis2[1])/np.pi,np.array(Phis2[0])/np.pi-np.array(Phis2[2])/np.pi,"brown"]
Zeta3 = [np.array(Phis3[0])/np.pi-np.array(Phis3[1])/np.pi,np.array(Phis3[0])/np.pi-np.array(Phis3[2])/np.pi,"green"]
#"""
plt.plot(R1, c = "blue")
plt.plot(R2, c="brown")
plt.show()
plt.plot(R3, c = "green")
plt.plot(R4, c = "red")
plt.show()
plt.plot(np.linspace(0,6000,60000),R1[:60000], c = "blue")
plt.plot(np.linspace(0,6000,60000),R2[:60000], c="brown")
plt.plot(np.linspace(0,6000,60000),R3[:60000], c = "green")
plt.plot(np.linspace(0,6000,60000),R4[:60000], c = "red")
plt.show()
"""
X = np.linspace(0,6000,60000)
plt.plot(X,Zeta3[0][:60000], c="green")
plt.plot(X,Zeta3[1][:60000], c="green", dashes = [5,2,2,2])
plt.show()
Zetas = [Zeta,Zeta2,Zeta3,Zeta4]
for Angle in Zetas:
	plt.plot(X,Angle[0][:60000], c = Angle[2])
	plt.plot(X,Angle[1][:60000], c = Angle[2], dashes=[5,2,2,2])
plt.show()
plt.plot(Zeta3[0][:60000],Zeta3[1][:60000], c = "green")
plt.show()
"""
