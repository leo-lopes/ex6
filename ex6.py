#exercicio 6 

def gaussian(A, w, x0, p,xs):
    u = ((A*w)/w)*np.exp(-(xs-x0)**2/(2*w**2) + (1j)*p*(xs-x0)) 
    u = u / np.sqrt((np.abs(u)**2).sum())
    return u
def gaussian2(A, w, x0, p,xs):
    u = (A*np.exp(-(xs-x0)**2/(2*w**2))) 
    return u
def sol(A, w, x0, p,xs,t):
    u = ((A*w)/np.sqrt(w**2 + 2*(1j)*t))*np.exp(-(xs-x0-2*p*t)**2/(2*w**2 + 4*(1j)*t) + (1j)*p*(xs-x0) - (1j)*(p**2)*t) 
    u = u / np.sqrt((np.abs(u)**2).sum())
    return u
def probability(state):
    prob = np.abs(state)**2
    return prob / prob.sum()
def potential (x):
    if x > 300 and x<400:
        V=np.sin(x/20.0*L)
    else:
        V=0
    return V
import numpy as npa
import matplotlib.pyplot as plt
import os, sys
import matplotlib

matplotlib.rc('font', size=18)
matplotlib.rc('font', family='Arial')
N = 501 #
dt = 5.e-2 
L = float(500)
nsteps = 50000 
dx = L/(N-1) 
nplot = 80 

alpha = ((1j)*dt)/(2*dx**2) 
A = np.zeros((N,N),dtype='complex')
B = np.zeros((N,N),dtype='complex')
for i in range(N):
    if i==0:
        A[i,:] = [1+2*alpha + (1j)*dt*potential(i*dx)/2.  if j==0 else -alpha if j==1 or j==N-1 else 0 for j in range(N)]
        B[i,:] = [1-2*alpha - (1j)*dt*potential(i*dx)/2. if j==0 else alpha if j==1 or j==N-1 else 0 for j in range(N)]
    elif i==N-1:
        A[i,:] = [-alpha if j==0 or j==N-2 else 1+2.*alpha + (1j)*dt*potential(i*dx)/2. if j==N-1 else 0 for j in range(N)]
        B[i,:] = [alpha if j==0 or j==N-2 else 1-2.*alpha - (1j)*dt*potential(i*dx)/2. if j==N-1 else 0 for j in range(N)]
    else:
        A[i,:] = [-alpha if j==i-1 or j==i+1 else 1+2.*alpha + (1j)*dt*potential(i*dx)/2. if j==i else 0 for j in range(N)]
        B[i,:] = [alpha if j==i-1 or j==i+1 else 1-2.*alpha - (1j)*dt*potential(i*dx)/2. if j==i else 0 for j in range(N)]
x = np.linspace(0, L, N)
t= np.arange(0,dt*nsteps,dt)
#initial condition
u = gaussian(1.0, L/50.0, L/5.,2.0, x) #evaluate right hand side at t=0
bb = B.dot(u[:]) 
sol=np.asarray([sol(1.0, L/50.0, L/5.,2.0, x,xx) for xx in t])
inv_A = np.linalg.inv(A)
c = 0
fig=plt.figure(figsize=(25.0,15.0))
plt.xlim(left=0,right=L)
plt.ylim(bottom=-1.1,top=1.1)
plt.tick_params(axis='x', labelsize=24)
plt.tick_params(axis='y', labelsize=24)
plt.title('(a) Simulação', size=24)
plt.axvspan(300., 400., facecolor='b', alpha=0.5)
plt.plot(x,u.real,linewidth=2)
plt.plot(x,u.imag,linewidth=2)
#plt.plot(x,pot,'r','--',linewidth=2)
plt.plot(x,probability(u),linewidth=2)
plt.legend(['Real','imag','prob'],prop={'size':10})
filename = 'foo' + str(c+1).zfill(3) + '.jpg';
plt.savefig(filename)
plt.clf()
for j in range(nsteps):
    #find solution inside domain
    u[:] = inv_A @ bb
    v=0.0
    for k in range(len(u)):
        v = v +abs(u[k])
    print(v)
    print(j)
    #update right hand side
    bb = B.dot(u[:]) 
    if(j%nplot==0): #plot results every nplot timesteps
        fig=plt.figure(figsize=(25.0,15.0))
        plt.xlim(left=0,right=L)
        plt.ylim(bottom=-1.1,top=1.1)
        plt.tick_params(axis='x', labelsize=24)
        plt.tick_params(axis='y', labelsize=24)
        plt.title('(a) Simulação', size=24)
        plt.axvspan(300., 400., facecolor='b', alpha=0.5)
        plt.plot(x,u.real,linewidth=2)
        plt.plot(x,u.imag,linewidth=2)
        plt.plot(x,probability(u),linewidth=2)
        #plt.plot(x,pot,'r','--',linewidth=2)
        plt.legend(['Real','imag','prob'],prop={'size':10})
        filename = 'foo' + str(c+1).zfill(3) + '.jpg';
        plt.savefig(filename)
        plt.clf()
        c += 1
os.system("ffmpeg -r 5 -y -i 'foo%03d.jpg' testando_novo_code.m4v")
os.system("rm -f *.jpg")
