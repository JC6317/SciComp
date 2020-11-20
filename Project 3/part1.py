"""Scientific Computation Project 3, part 1
Your CID here:01063446
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy as sp

def hfield(r,th,h,levels=50):
    """Displays height field stored in 2D array, h,
    using polar grid data stored in 1D arrays r and th.
    Modify as needed.
    """
    thg,rg = np.meshgrid(th,r)
    xg = rg*np.cos(thg)
    yg = rg*np.sin(thg)
    plt.figure()
    plt.contourf(xg,yg,h,levels)
    plt.axis('equal')
    #change colours
    return None

def repair1(R,p,l=1.0,niter=10,inputs=()):
    """
    Question 1.1: Repair corrupted data stored in input
    array, R.
    Input:
        R: 2-D data array (should be loaded from data1.npy)
        p: dimension parameter
        l: l2-regularization parameter
        niter: maximum number of iterations during optimization
        inputs: can be used to provide other input as needed
    Output:
        A,B: a x p and p x b numpy arrays set during optimization
    """
    #problem setup
    R0 = R.copy()
    a,b = R.shape
    iK,jK = np.where(R0 != -1000) #indices for valid data
    aK,bK = np.where(R0 == -1000) #indices for missing data

    S = set()
    for i,j in zip(iK,jK):
            S.add((i,j))

    #Set initial A,B
    A = np.ones((a,p))
    B = np.ones((p,b))

    #Create lists of indices used during optimization
    mlist = [[] for i in range(a)]
    nlist = [[] for j in range(b)]

    for i,j in zip(iK,jK):
        mlist[i].append(j)
        nlist[j].append(i)

    dA = np.zeros(niter)
    dB = np.zeros(niter)

    for z in range(niter):
        Aold = A.copy()
        Bold = B.copy()

        #Loop through elements of A and B in different
        #order each optimization step
        for m in np.random.permutation(a):
            for n in np.random.permutation(b):
                if n < p: #Update A[m,n]
                    Bfac = 0.0
                    Asum = 0

                    for j in mlist[m]:
                        Bfac += B[n,j]**2
                        Rsum = 0
                        for k in range(p):
                            if k != n: Rsum += A[m,k]*B[k,j]
                        Asum += (R[m,j] - Rsum)*B[n,j]

                    A[m,n] = Asum/(Bfac+l) #New A[m,n]
                if m<p:
                    #Add code here to update B[m,n]
                    B[m,n]=None #modify
        dA[z] = np.sum(np.abs(A-Aold))
        dB[z] = np.sum(np.abs(B-Bold))
        if z%10==0: print("z,dA,dB=",z,dA[z],dB[z])


    return A,B


def repair2(R,p,l=1.0,niter=10,inputs=()):
    """
    Question 1.1: Repair corrupted data stored in input
    array, R. Efficient and complete version of repair1.
    Input:
        R: 2-D data array (should be loaded from data1.npy)
        p: dimension parameter, chosen rank?
        l: l2-regularization parameter
        niter: maximum number of iterations during optimization
        inputs: can be used to provide other input as needed
    Output:
        A,B: a x p and p x b numpy arrays set during optimization
    """
    #problem setup
    R0 = R.copy()
    a,b = R.shape
    iK,jK = np.where(R0 != -1000) #indices for valid data
    aK,bK = np.where(R0 == -1000) #indices for missing data

    S = set()
    for i,j in zip(iK,jK):
            S.add((i,j))

    #Set initial A,B
    A = np.ones((a,p))
    B = np.ones((p,b))
    
    dA = np.zeros(niter)
    dB = np.zeros(niter)
    
    mdict = {k:[] for k in range(300)}#create dictionary to find i and j vectors later
    ndict = {k:[] for k in range(300)}
    for i,j in zip(iK,jK):
        mdict[i].append(j)
        ndict[j].append(i)

    #iterations defined from input
    for z in range(niter):
        print('iteration',z,'out of',niter)
        Aold = A.copy()
        Bold = B.copy()

        #Loop through elements of A and B in random order
        for m in np.random.permutation(a): #in iK?
            #j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)
            j=mdict[m]
            for n in np.random.permutation(b):
                #j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)

                if n<p:
    
                    #j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)
                    #np.delete needs to specify 0 or 1 for row or column delete, otherwise array is unpacked into list
                    Rmj_less_sum = R[m,j]- np.dot(np.delete(A[m,:],n) , np.delete(B[:,j],n,0)) #inner product of Amk,Bkj, k is not n
                    entry = np.dot(Rmj_less_sum,B[n,j]) #j is 1d, so B[n,j] is a vector, so this is inner product
                    entry = entry/(l + np.dot(B[n,j],B[n,j]))#divided by sum on LHS
                    A[m,n] = entry
                
                if m<p:
                    #i=iK[np.where(jK==n)] #i needs to be vector of entries where i,n is valid entry (i,n) in (iK,jK)
                    i=ndict[n]
                    Rin_less_sum = R[i,n]- np.dot(np.delete(A[i,:],m,1) , np.delete(B[:,n],m)) #inner product of Aik,Bkn, k is not m
                    entry2 = np.dot(Rin_less_sum,A[i,m]) #i is 1d, so A[i,m] is a vector, so this is inner product
                    entry2 = entry2/(l + np.dot(A[i,m],A[i,m]))#divided by sum on LHS
                    B[m,n] = entry2
                    
        #checking the change in updated A and B
        dA[z] = np.sum(np.abs(A-Aold))
        dB[z] = np.sum(np.abs(B-Bold))
        if z%10==0: print("z,dA,dB=",z,dA[z],dB[z])
                    
    return A,B


def outwave(r0):
    """
    Question 1.2i)
    Calculate outgoing wave solution at r=r0
    See code/comments below for futher details
        Input: r0, location at which to compute solution
        Output: B, wave equation solution at r=r0

    """
    A = np.load('data2.npy')
    r = np.load('r.npy')
    th = np.load('theta.npy')

    Nr,Ntheta,Nt = A.shape
    B = np.zeros((Ntheta,Nt))

    return B

def analyze1():
    """
    Question 1.2ii)
    Add input/output as needed

    """
    
    #array is 2pi grid, so pi/4 is 1/8 into grid, theta has 289 values
    pi14 = 288//8 #this is the index where th[pi14] = pi/4
    pi34 = 288*3//8 #index for 3pi/4
    pi54 = 288*5//8 #index for 5pi/4
    
    
    wave_1 = data3[:,pi14,:] #slice at pi/4 
    wave_3 = data3[:,pi34,:] #slice at 3pi/4
    wave_5 = data3[:,pi54,:] #slice at 5pi/4
    
    #plt.plot(wave_1)    
    #plt.plot(wave_3)
    #plt.plot(wave_5)


    
    def fix_t_plot(t):#function to plot waves (height against radial distance) at fixed time 
        plt.plot(wave_1[:,t], label='pi/4 wave') #wave_1 at fixed time
        plt.plot(wave_3[:,t], label='3pi/4 wave') #wave_3 at fixed time
        plt.plot(wave_5[:,t], label='5pi/4 wave') #wave_5 at fixed time
        plt.xlabel('ith radial point (from 1 to 5)')
        plt.ylabel('wave amplitude')
        plt.legend(loc="lower right")
        plt.title('Waves for different theta values at a fixed time: %i th time point' % t)
        plt.show()
        means = np.mean(wave_1[:,t]), np.mean(wave_3[:,t]), np.mean(wave_5[:,t])
        varis = np.var(wave_1[:,t]), np.var(wave_3[:,t]), np.var(wave_5[:,t])

        print('pi/4 mean, 3pi/4 mean, 5pi/4 mean:', means, 'pi/4 var, 3pi/4 var, 5pi/4 var', varis)
    

    def fft_tplot(wave,t,Nt=300):    #fourier transform coefficients for a fixed time
        w = wave[:,t] #extract wave at time
        c = np.fft.fft(w) #fft coefficients
        c = np.fft.fftshift(c)/(Nt/2) #rearrange so ordered
        n = np.arange(-Nt/2,Nt/2)
        plt.figure()
        plt.plot(n,np.abs(c), 'x')
        plt.xlim(0,150) #this limit halves the graphs since the functions are real valued and we do not need duplicate coefficients
        plt.xlabel('mode number, n')
        plt.ylabel('$|c_n|$')
        plt.title('Fourier Transform coefficients for wave5 at time %i' %t)
    
    

    fix_t_plot(0)
    plt.figure()
    
    fft_tplot(wave_1,0)
    plt.figure()

    fft_tplot(wave_3,0)
    plt.figure()

    fft_tplot(wave_5,0)
    plt.figure()
    
    fix_t_plot(60)
    plt.figure()

    fft_tplot(wave_1,60)
    plt.figure()

    fft_tplot(wave_3,60)
    plt.figure()

    fft_tplot(wave_5,60)
    
    return None #modify as needed
    

def reduce(H,inputs=12):
    """
    Question 1.3: Construct one or more arrays from H
    that can be used by reconstruct
    Input:
        H: 3-D data array
        inputs: can be used to provide other input as needed
    Output:
        arrays: a tuple containing the arrays produced from H
    """
    """
    Method: convert 3d array to 2d array. Apply PCA technique using SVD 
    """
    x,y,z = H.shape #get dimension of H where H is 3d
    
    #H_2d = H.transpose(2,0,1).reshape(-1,H.shape[1]) #unpack H into 2d by z dimension
    H_2d = np.reshape(H,(x*z,y)) #unpack into 2d default way
    
    M,N = np.shape(H_2d) #get dimensions of this 2d data array
    #print("M,N=",M,N)
    H_2d_mr = H_2d - np.outer(np.ones((M,1)),H_2d.mean(axis=0)) #H with mean removed
    #print(np.mean(H_2d_mr,axis=0))
    
    H_BSR_T = sp.sparse.bsr_matrix(H_2d_mr.T) #sparse conversion otherwise memory error (transposed)
    U,S,VT = sp.sparse.linalg.svds(H_BSR_T,k=inputs) #svd on transposed sparse 2d matrix with mean removed
    G = np.dot(U.T,H_2d_mr.T) #G= Ut * A, where A is the transpose of mean removed data 
         
    arrays = (U,G,np.array([x,y,z])) #tuple of arrays
    return arrays

#Reshape2 = R2.transpose(2,0,1).reshape(-1,R2.shape[1])

def reconstruct(arrays,inputs=()):
    """
    Question 1.3: Generate matrix with same shape as H (see reduce above)
    that has some meaningful correspondence to H
    Input:
        arrays: tuple generated by reduce
        inputs: can be used to provide other input as needed
    Output:
        Hnew: a numpy array with the same shape as H
    """
    U,G,shape = arrays #unpack arrays
    Hprime_2d = np.dot(U,G) #UT is the transformation vector and is orthogonal, so U is inverse
    Hprime= Hprime_2d.T #tranpose before reshaping
    Hnew= Hprime.reshape(shape[0],shape[1],shape[2]) #reshape back to original dimensions

    return Hnew


if __name__=='__main__':
    x=None
    #Add code here to call functions above and
    #generate figures you are submitting
    #import os
    #os.chdir('C:\\Users\\Jon\\Documents')
    #os.chdir('C:\\Users\\JC\\Documents')
    
    data1 = np.load('data1.npy')
    r=np.load('r.npy')
    th=np.load('theta.npy')
    
    #q1
    hfield(r,th,data1) #plot data1
    #np.linalg.matrix_rank(data1) #rank of data1 is what we will set p as (36)
    np.random.seed(112)
    A_prime,B_prime = repair2(data1,36,niter=101,l=1)
    data1_rest = A_prime @ B_prime
    hfield(r,th,data1_rest)
    
    data3 = np.load('data3.npy')
    
    plt.figure()
    analyze1()
    
    plt.figure()
    arrays = reduce(data3)
    data3new = reconstruct(arrays)
    hfield(r,th,data3[:,:,0])
    hfield(r,th,data3new[:,:,0])

