""" CID:01063446:
    Template code for part 1, contains 4 functions:
    newSort, merge: codes for part 1.1
    time_newSort: to be completed for part 1.1
    findTrough: to be completed for part 1.2
"""
import numpy as np
import timeit
import time
import matplotlib.pyplot as plt

def newSort(X,k=0):
    """Given an unsorted list of integers, X,
        sort list and return sorted list
    """
    n = len(X)
    if n==1:
        return X
    elif n<=k:
        for i in range(n-1):
            ind_min = i
            for j in range(i+1,n):
                if X[j]<X[ind_min]:
                    ind_min = j
            X[i],X[ind_min] = X[ind_min],X[i]
        return X
    else:
        L = newSort(X[:n//2],k)
        R = newSort(X[n//2:],k)
        return merge(L,R)
"""test correctness"""
def newSort_test(length,k=0): #length of list, threshold for merge
    #length = np.random.randint(0,1000)
    X=np.random.randint(-2*length,2*length,length) #
    X_sorted = newSort(X,k)
    check = True #assume check is true
    for i in range(X.shape[0]-1): #-1 because we can not compare the last element of the list with anything on the right
        if X_sorted[i]>X_sorted[i+1]:
            check = False #check changes to false if an integer is greater than one it should less equal to
            break
    return check

"""batch test correctness"""
def newSort_batchtest(tests): #tests is number of tests to run
    results = np.zeros(tests)
    for t in range(tests):
        length = np.random.randint(1,100) #generate random length from 1-1000
        k = np.random.randint(1,length) #generate random k
        check = newSort_test(length,k) #test whether sort has worked correctly
        results[t] = check
    return str(np.sum(results))+" successes out of "+str(tests) #sum counts no. of trues

"""test sorting time for different lengths"""
def newSort_time_length(l_start,l_end,interv,k=0,average=1): #starting length, ending length, no. of intervals
    results = np.zeros((interv+1,1+average)) #preallocate empty array to store times, first column stores length, subsequent columns stores times from a loop
    for a in range(1,average+1,1): #each loop is a new run of varying lengths    
        for i in range(interv+1):
            length = l_start + i*(l_end-l_start)//interv
            X=np.random.randint(1,2*length,length) #low, high, length of list
            t1 = time.time()
            newSort(X,k) #function runs here
            t2 = time.time() #time
            
            results[i,0]=length #storing the length in first column
            results[i,a]=t2-t1 #storing the time is a column
            
    plt.plot(results[:,0],np.mean(results[:,1:],axis=1),label='k=%1.0f'%(k))
    plt.xlabel("list length")
    plt.ylabel("time (seconds)")
    plt.legend(loc='upper left')
    #plt.savefig('figlength2000_k%1.0f.png'%(k)) #saving figures
    #plt.savefig("figlength2000_maxk_noaverage")
    return results
 
"""test sorting time for different k"""
def newSort_time_threshold(k_start,k_end,interv,length=10000): #starting length, ending length, no. of intervals
    results = np.zeros((interv+1,2)) #preallocate empty array to save times
    X=np.random.randint(1,2*length,length) #list goes outside loop since we want to vary k for the same list
    for i in range(interv+1):
        k = k_start + i*(k_end-k_start)//interv
        #X=np.random.randint(1,2*length,length) #low, high, length of list
        t1 = time.time()
        newSort(X,k)
        t2 = time.time()
        
        results[i,0]=k
        results[i,1]=t2-t1
    
    plt.plot(results[:,0],results[:,1])
    plt.xlabel("K value")
    plt.ylabel("time (seconds)")
    #plt.savefig("fig_threshold_k_end%1.0f.png"%(k_end))
    return results

def merge(L,R):
    """Merge 2 sorted lists provided as input
    into a single sorted list
    """
    M = [] #Merged list, initially empty
    indL,indR = 0,0 #start indices
    nL,nR = len(L),len(R)

    #Add one element to M per iteration until an entire sublist
    #has been added
    for i in range(nL+nR):
        if L[indL]<R[indR]:
            M.append(L[indL])
            indL = indL + 1
            if indL>=nL:
                M.extend(R[indR:])
                break
        else:
            M.append(R[indR])
            indR = indR + 1
            if indR>=nR:
                M.extend(L[indL:])
                break
    return M


"""test sorting time for either length or k"""
def time_newSort(variety,start,end,interv,k=0,length=10000,average=1):
    """Analyze performance of newSort
    Use variables inputs and outputs if/as needed
    """
    if variety==0: #variety==0 calls the function that varies length for fixed k
        outputs = newSort_time_length(start, end, interv,k,average) #uses the inputs from above
        
    elif variety==1: #variety==1 calls the function that varies k for fixed length
        outputs = newSort_time_threshold(start, end, interv,length)
        
    else:
        print("set variety to either 0 or 1, to vary length or k-value respectively")     
    
    return outputs


def findTrough(L,start=0,end=-1): #input a list L, start and end are not set
    """Find and return a location of a trough in L
    """
    if end ==-1: #will only trigger on first run
        end = len(L)-1 #sets length to starting length
        
    if start>end: #triggers when trough is not found
        return -1 #output when no trough is found

    
    mid = (start + end)//2 #this sets the index to the midpoint of the list, integer division is used so both odd/even lengths work
    
    if start == mid or end==mid: #this condition is needed in case like L=[3,2,1] or L=[1,2,3,4,5,0]
        return mid
    
    if L[mid]<= L[mid-1] and L[mid]<=L[mid+1]: #this is the condition for a trough
        return mid #function ends here is trough is found
    elif L[mid]>L[mid-1]: #searching the left half
        end = mid -1 #set a new end/right point
        return findTrough(L,start,end) #recursive call
    else: #L[mid]>L[mid+1]
        start = mid +1
        return findTrough(L,start,end)
        
        


if __name__=='__main__':
    #inputs=None
    #outputs=time_newSort(inputs) 
    
    newSort_batchtest(100) #testing the accuracy of newSort for 100 random lists
    
    #newSort_time_length(1,2001,200,0,1)
    time_newSort(0,1,2001,200,k=20)  #non averaged run for k=0 varying length
    
    for i in [0,5,10,20,40]: #looking for the optimum k value, plotting a series of graphs with k=0,5,10,20,40 and varying length for each
        np.random.seed(1)
        time_newSort(0,1,2001,200,k=i,average=10)
    
    time_newSort(0, 1, 2001,200,k=2000,average=10) #setting k=2000 and varying length

    time_newSort(1,0,350,70) #varying k from 0 to 350 for fixed length list
    
    for i in range(5): #varying k from 0-20 for fixed length list and plotting multiple graphs
        time_newSort(1,0,20,20)