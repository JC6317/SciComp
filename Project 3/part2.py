"""Scientific Computation Project 3, part 2
Your CID here:01063446
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy as sp
import time as ti

def microbes(phi,kappa,mu,L = 1024,Nx=1024,Nt=1201,T=600,display=False):
    """
    Question 2.2
    Simulate microbe competition model

    Input:
    phi,kappa,mu: model parameters
    Nx: Number of grid points in x
    Nt: Number of time steps
    T: Timespan for simulation is [0,T]
    Display: Function creates contour plot of f when true

    Output:
    f,g: Nt x Nx arrays containing solution
    """

    #generate grid
    L = 1024
    x = np.linspace(0,L,Nx)
    dx = x[1]-x[0]
    dx2inv = 1/dx**2

    def RHS(y,t,k,r,phi,dx2inv):
        #RHS of model equations used by odeint

        n = y.size//2

        f = y[:n]
        g = y[n:]

        #Compute 2nd derivatives
        d2f = (f[2:]-2*f[1:-1]+f[:-2])*dx2inv
        d2g = (g[2:]-2*g[1:-1]+g[:-2])*dx2inv

        #Construct RHS
        R = f/(f+phi)
        dfdt = d2f + f[1:-1]*(1-f[1:-1])- R[1:-1]*g[1:-1]
        dgdt = d2g - r*k*g[1:-1] + k*R[1:-1]*g[1:-1]
        dy = np.zeros(2*n)
        dy[1:n-1] = dfdt
        dy[n+1:-1] = dgdt

        #Enforce boundary conditions
        a1,a2 = -4/3,-1/3
        dy[0] = a1*dy[1]+a2*dy[2]
        dy[n-1] = a1*dy[n-2]+a2*dy[n-3]
        dy[n] = a1*dy[n+1]+a2*dy[n+2]
        dy[-1] = a1*dy[-2]+a2*dy[-3]

        return dy


    #Steady states
    rho = mu/kappa
    F = rho*phi/(1-rho)
    G = (1-F)*(F+phi)
    y0 = np.zeros(2*Nx) #initialize signal
    y0[:Nx] = F
    y0[Nx:] = G + 0.01*np.cos(10*np.pi/L*x) + 0.01*np.cos(20*np.pi/L*x)

    t = np.linspace(0,T,Nt)

    #compute solution
    print("running simulation...")
    y = odeint(RHS,y0,t,args=(kappa,rho,phi,dx2inv),rtol=1e-6,atol=1e-6)
    f = y[:,:Nx]
    g = y[:,Nx:]
    print("finished simulation")
    if display:
        plt.figure()
        plt.contour(x,t,f)
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('Contours of f')


    return f,g


def newdiff(f,h):
    """
    Question 2.1 i)
    Input:
        f: array whose 2nd derivative will be computed
        h: grid spacing
    Output:
        d2f: second derivative of f computed with compact fd scheme
    """
    N = f.shape[0]
    hfac = 1/(h*h) #change
    
    #Coefficients for compact fd scheme
    alpha = 9/38
    a = (696-1191*alpha)/428
    b = (2454*alpha-294)/535
    c = (1179*alpha-344)/2140
    
    """
    our finite difference scheme has 3 equations
    f''_0 + 10f''_1 =...upper boundary UB
    alphaf''_i-1 + f''_i + alphaf''_i+1 =... which is used for the MEAT of this scheme
    f''_N-1 + 10f''_N-2 =...lower boundary LB
    
    we will calculate a single UB and LB value
    and calculate a 1d array of MEAT values
    to create an array b = [UB,MEAT,LB]
    Now each of these b entries corresponds to a a lin-combo of 2nd derivs
    so we create a banded matrix A, where the the jth column corresponds to coefficients 2nd deriv at jth point in grid
    and the ith row represents the the choice of i for the equations above
    x is a vector of actual values of the derivative 
    we use a banded solver to solve this system of equations to return x, our d2f
    """
    u0,u1,u2,u3,u4 = 145/12,76/3,29/2,4/3,1/12 #set coefficients for UB equation for neatness
    ln1,ln2,ln3,ln4,ln5 = 145/12,76/3,29/2,4/3,1/12 #coefficients for LB equations (ln1:L_N-1 etc)
    
    UB = hfac*(u0*f[0] - u1*f[1] + u2*f[2] - u3*f[3] + u4*f[4]) #upper boundary and first value of b vector
    LB = hfac*(ln1*f[N-1] - ln2*f[N-2] + ln3*f[N-3] - ln4*f[N-4] +ln5*f[N-5]) #lower boundary and last value of b vec
    
    #creating indexes, later we will use .take with mode='wrap' so that we can call elements of the element cyclically 
    in3 = range(1-3,N-1-3) #note that f is periodic
    i3 = range(1+3,(N-1+3))
    in2 = range(1-2,N-1-2)
    i2 = range(1+2,(N-1+2))
    in1 = range(1-1,N-1-1) #the -1 and =1 could be done without wrap (potentially faster), but we will do all the slicing in the same way for consistency
    i1 = range(1+1,(N-1+1))    #maybe create a function later for neatness?
    meat =  (c/9)*(f.take(in3,mode='wrap')- 2*f[1:N-1] + f.take(i3,mode='wrap'))\#wrap allows us to to cycle from one end of array to the other
        + (b/4)*(f.take(in2,mode='wrap')- 2*f[1:N-1] + f.take(i2,mode='wrap'))\
            + a *(f.take(in1,mode='wrap')- 2*f[1:N-1] + f.take(i1,mode='wrap')) #this is the meat of the b vector
    
    bvec = np.hstack((UB,hfac*meat,LB)) #adding the values from boundary conditions to either end of meat

    diag_0 = np.ones(N) #create diagonals as rows 
    diag_1 = np.hstack((10,alpha*np.ones(N-2)))
    diag_n1= np.hstack((alpha*np.ones(N-2),10)) 
    AB = np.array([np.hstack((0,diag_1)),diag_0,np.hstack((diag_n1,0))]) #stack rows in matrix diagonal ordered form
    (l,u) = (1,1)
    #solve banded for Ax = b, x=d2f
    d2f = sp.linalg.solve_banded((l,u),AB,bvec) 
    
    return d2f #modify as needed

def analyzefd():
    """
    Question 2.1 ii)
    Add input/output as needed

    """
    N=1000#mumber of times to take
    
    results=np.zeros(N)#preallocate space for results
    Nlengths = np.linspace(10,90000,N) #create grid of points for function
    for i in range(N):#
        N=Nlengths[i]
        x = np.linspace(0,2*np.pi, int(N))
        fsin = np.sin(x)
        start = ti.time() #time start
        d2fsin = newdiff(fsin,2*np.pi/N) #apply newdiff
        end = ti.time() #time end
        results[i] = end-start
    
    plt.plot(Nlengths,results)
    plt.title('newdiff runtime against f length')
    plt.xlabel('length of f array')
    plt.ylabel('runtime (seconds)')
    
    
    n=np.linspace(0,2*np.pi,100)
    d2f = newdiff(np.sin(n),2*np.pi/100) #calculated deriv
    true_d2f = -np.sin(n) #second erivative of sin is -sin
    
    
    plt.figure()
    error = d2f -true_d2f #difference between calculated and theoretical
    plt.plot(range(100),error)
    plt.title('difference between true second derivative and newdiff output')
    plt.xlabel('ith point in array')
    plt.ylabel('error')
    
    return None #modify as needed


def dynamics(f,m,elo=0.1,ehi=20,points=1000):
    """
    Question 2.2
    Add input/output as needed
    """
    
    ft = f[200:,:] #discard transient
    D = sp.spatial.distance.pdist(ft)
    epsi = np.linspace(elo,ehi,points) #log not lin!
    
    C=np.zeros_like(epsi)
    for i in range(len(epsi)):
        Di = D[D<epsi[i]] #discard D larger than epsilon
        C[i] =Di.size #new correlation sum unscaled
    
    C = C*m*(m-1)/2 #scaling
    plt.plot(epsi,C)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(elo,ehi)
    plt.xlabel('epsilon')
    plt.ylabel('C')
    plt.title('correlation sum plot')
    

    return np.polyfit(np.log(epsi),np.log(C),1)
    
    
if __name__=='__main__':
    x=np.linspace(0,2*np.pi,1000)
    f=np.sin(x)
    d2f =newdiff(f,2*np.pi/1000)
    plt.plot(d2f[4:1000-4],label='numerical d2f(sin)')
    plt.plot(np.sin(x[4:996]),label='sine function')
    plt.title('plot of sin function against 2nd derivative, with the outside points excluded')
    plt.legend()
    plt.figure()
    analyzefd()
    
    plt.figure()
    f2,g2 = microbes(0.3,2,0.8)
    dynamics(f2,1024,0.1,20,1000)
    plt.figure()
    dynamics(f2,1024,2,6,1000)
    
    plt.figure()
    f17,g17 = microbes(0.3,1.7,0.4*1.7)
    dynamics(f17,1024,2,6,1000)

    plt.figure()    
    f15,g15 = microbes(0.3,1.5,0.4*1.5)
    dynamics(f15,1024,2,6,1000)



