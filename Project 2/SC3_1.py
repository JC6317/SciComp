"""Scientific Computation Project 3, part 1
Your CID here
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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
    return None

r = np.linspace(1,5,300) #300 is R rows
th = np.linspace(0,2*np.pi,289)

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
        for m in np.random.permutation(a): #in iK?
            for n in np.random.permutation(b):

                j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)
    
                Rmj_less_sum = R[m,j]- np.dot(np.delete(A[m,:],n) , np.delete(B[:,j],n)) #inner product of Amk,Bkj, k is not n
                entry = np.dot(Rmj_less_sum,B[n,j]) #j is 1d, so B[n,j] is a vector, so this is inner product
                entry = entry/(np.dot(B[n,j],B[n,j]))#divided by sum on LHS
                A[m,n] = entry
                
                i=iK[np.where(jK==n)] #i needs to be vector of entries where i,n is valid entry (i,n) in (iK,jK)

                Rin_less_sum = R[i,n]- np.dot(np.delete(A[i,:],m) , np.delete(B[:,n],m)) #inner product of Aik,Bkn, k is not m
                entry2 = np.dot(Rin_less_sum,A[i,m]) #i is 1d, so A[i,m] is a vector, so this is inner product
                entry2 = entry/(np.dot(A[i,m],A[i,m]))#divided by sum on LHS
                B[m,n] = entry
                    
                
                
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

    #iterations defined from input
    for z in range(niter):
        print(z)
        Aold = A.copy()
        Bold = B.copy()

        #Loop through elements of A and B in different
        #order each optimization step
        for m in np.random.permutation(a): #in iK?
            j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)
            for n in np.random.permutation(b):
                #print('n',n)
                #j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)

                if n<p:
    
                    #j=jK[np.where(iK==m)] #j needs to be vector of entries where m,j is valid entry (m,j) in (iK,jK)
                    #np.delete needs to specify 0 or 1 for row or column delete, otherwise array is unpacked into list
                    Rmj_less_sum = R[m,j]- np.dot(np.delete(A[m,:],n) , np.delete(B[:,j],n,0)) #inner product of Amk,Bkj, k is not n
                    entry = np.dot(Rmj_less_sum,B[n,j]) #j is 1d, so B[n,j] is a vector, so this is inner product
                    entry = entry/(np.dot(B[n,j],B[n,j]))#divided by sum on LHS
                    A[m,n] = entry
                
                #print('m',m)
                if m<p:
                    i=iK[np.where(jK==n)] #i needs to be vector of entries where i,n is valid entry (i,n) in (iK,jK)
    
                    Rin_less_sum = R[i,n]- np.dot(np.delete(A[i,:],m,1) , np.delete(B[:,n],m)) #inner product of Aik,Bkn, k is not m
                    entry2 = np.dot(Rin_less_sum,A[i,m]) #i is 1d, so A[i,m] is a vector, so this is inner product
                    entry2 = entry2/(np.dot(A[i,m],A[i,m]))#divided by sum on LHS
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

    return None #modify as needed




def reduce(H,inputs=()):
    """
    Question 1.3: Construct one or more arrays from H
    that can be used by reconstruct
    Input:
        H: 3-D data array
        inputs: can be used to provide other input as needed
    Output:
        arrays: a tuple containing the arrays produced from H
    """

    #Add code here
    x,y,z = H.shape #get dimension of H where H is 3d
    
    #H_2d = H.transpose(2,0,1).reshape(-1,H.shape[1]) #unpack H into 2d by z dimension
    H_2d = np.reshape(H,(x*z,y)) #unpack into 2d default way
    
    M,N = np.shape(H_2d) #get dimensions of this 2d data array
    print("M,N=",M,N)
    H_2d_mr = H_2d - np.outer(np.ones((M,1)),H_2d.mean(axis=0)) #H with mean removed
    print(np.mean(H_2d,axis=0))
    U,S,VT = np.linalg.svd(H_2d_mr.T) #svd
    
    G = np.dot(U.T,H_2d_mr.T) #G= Ut * A, where A is teh transpose of mean removed data 
    #first 3 rows of are dimension reduced
    #check svd to see howm any rows
    #return first few rows of G 
    #AND return U or V, possibly sparse
    
        
    arrays = () #modify as needed
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
    
    #Hnew = np.dot(V.T,PC)
    #Vt is matrix of evec from sample covariance of X
    #PC is output from above
    #need to add mean somewhere
    
    

    #Add code here

    Hnew = None #modify

    return Hnew


if __name__=='__main__':
    x=None
    #Add code here to call functions above and
    #generate figures you are submitting
    import os
    os.chdir('C:\\Users\\JC\\Documents')
    np.load('data1.npy')
