"""Scientific Computation Project 2, part 2
Your CID here:01063446
"""
import numpy as np
import networkx as nx
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import sparse,linalg


def nextboi(nodi,deci,G):
    """
    takes input node and decimal value, which should be between 0 and 1.
    finds neighbours of nodi, and total no. of neighbours.
    using the number of neighbours we can partition [0,1] into this many intervals
    np.ceil will find out which interval deci slots into and then returns this indexed neghibour
    NOTE!: This is to be run inside of rwgraph since rwgraph has an input toe define G
    edit:above updated
    """
    nb = list(G.neighbors(int(nodi))) #get neighbours
    index = int(np.ceil(deci*len(nb))) #int so we can index, ceil to find numbered interval
    #list comp to round?
    return nb[index-1]
    

def rwgraph(G,i0=0,M=100,Nt=100):
    """ Question 2.1
    Simulate M Nt-step random walks on input graph, G, with all
    walkers starting at node i0
    Input:
        G: An undirected, unweighted NetworkX graph
        i0: intial node for all walks
        M: Number of walks
        Nt: Number of steps per walk
    Output: X: M x Nt+1 array containing the simulated trajectories
    """
    X = np.zeros((M,Nt+1)) #zeros matrix to fill in with trajectories
    X[:,0]=i0 #set first column to be initial locaitons for each 

    #first step
    i0_nb = list(G.neighbors(i0)) #all neighbours of i0
    X[:,1] = np.random.choice(i0_nb,M) #first step can be easily calculated since all walkers start on i0
    
    #other steps
    X[:,2:] = np.random.uniform(0,1,(M,Nt-1)) #random numbers in [0,1] that will choose where each walker goes
    
    for c_i in range(2,Nt+1,1): #for each column index aka time step
        for r_i in range(M):
            X[r_i,c_i] = nextboi(X[r_i,c_i-1],X[r_i,c_i],G) #vectorised function that calculates next step for all M walkers
        
    return X

#np.random.uniform(0,1,(M,Nt-1))
#[(k,v) for k,v in deg_dict.items()][0]
#degrees = [val for (node, val) in G.degree()]
    #sorted(BAG_0.degree, key=lambda tup:tup[1])
"""
BAG_0 = nx.barabasi_albert_graph(2000,4,seed=0)
deg_arr2 = np.array(BAG_0.degree)
deg_dict = dict(BAG_0.degree)    
deg_dict_sort = {k: v for k, v in sorted(deg_dict.items(), key=lambda item: item[1])}
deg_dict_sort2 = {k: v for k, v in sorted(deg_dict.items(), key=lambda item: item[1],reverse=True)}
node0 = list(deg_dict_sort2)[0]

rwgraph(BAG_0,i0=node0, M=100,Nt=100)
#do for large M and Nt
#multiple runs measure average positions?
#plot line graph to show tending to limtiing distribution
#need large m because no of nodes

for a in range(averages):
    print(a)
    #rwgraph
    X = rwgraph(BAG_0,i0=node0,M=10000,Nt=200)
    #extract last column
    finalpos = X[:,X.shape[1]-1]
    #finalpos = X[:,Nt]
    uniqs,occurs = np.unique(finalpos,return_counts=True)
    #hist1 = [(node,deg) for node,deg in deg_dict.items()]
    plt.plot(hist2[:,0],hist2[:,1])
    hist2 = np.column_stack((uniqs,occurs))
    plt.plot(deg_arr2[:,0],deg_arr2[:,1])
    #pl.plot(hist[],hist[])
    #plot for 100 space steps and show convergence
    #increasing m
    
    #Q = sparse.diags(deg_arr2[:,1],format="csr")
    
    A = nx.adjacency_matrix(BAG_0)
    Q = sparse.diags(A.toarray().sum(axis=1),format="csr")
    
    L = Q - A
    e_val,e_vec = sparse.linalg.eigs(L.toarray(),k=2000)
    e_val,e_vec = np.linalg.eigh(L.toarray())
    
    
    invdeg = [(1/k) for k in deg_arr2[:,1]]
    L_s = invdeg * L
    e_val_s,e_vec_s = sparse.linalg.eigs(L_st.toarray(),k=2000)

    L_st= L_s.transpose
    e_val_st,e_vec_st = sparse.linalg.eigs(L_st.toarray(),k=2000)

"""   
    

def rwgraph_analyze1(walkers,time,deg_plot=True):
    """Analyze simulated random walks on
    Barabasi-Albert graphs.
    Modify input and output as needed.
    """
    BAG_0 = nx.barabasi_albert_graph(2000,4,seed=0) #fixed seed0
    X = rwgraph(BAG_0,i0=node0,M=walkers,Nt=time)
    finalpos = X[:,time]
    uniqs,occurs = np.unique(finalpos,return_counts=True)
    hist2 = np.column_stack((uniqs,occurs))
    plt.plot(hist2[:,0],hist2[:,1])
    plt.xlabel('node')
    plt.ylabel('frequency')
    plt.title('Random walker node frequencies at t='+str(time)+',M='+str(walkers))
    
    if deg_plot: #do you want to plot the degrees
        deg_arr2 = np.array(BAG_0.degree)
        plt.figure(1)
        plt.plot(deg_arr2[:,0],deg_arr2[:,1])
        plt.xlabel('node')
        plt.ylabel('degree')
        plt.title('Node degrees for BA m=2000,n=4,seed=0')

    return None #modify as needed


def rwgraph_analyze2(input=(None)):
    """Analyze similarities and differences
    between simulated random walks and linear diffusion on
    Barabasi-Albert graphs.
    Modify input and output as needed.
    """
    BAG_0 = nx.barabasi_albert_graph(2000,4,seed=0)

    A = nx.adjacency_matrix(BAG_0)
    Q = sparse.diags(A.toarray().sum(axis=1),format="csr")
    Q0 = A.toarray().sum(axis=1)
    node0=np.argmax(Q0)

    L = Q - A
    L = L.toarray()
    alpha,t,t0=(0.1,1000,0.1)
    it_0 = np.zeros(2000)
    it_0[node0] = t0
    
    res = []
    
    for t in np.linspace(0,5,1000):#sparse efficiency warnings
        res.append(linalg.expm(-L*t)@it_0)
        
    plot(res)
        
    
    
    it_0 = sparse.csr_matrix(it_0) #i of t0
    i_t = linalg.expm(-alpha * L * t)@ it_0 #too slow to run
    #suppose to take matrix exponential of -alpha *L *t with large t
    
    

    Qinv = sparse.linalg.inv(Q) #change to spsolve and solve for identity
    
    L_s = Qinv * L
    e_val_s,e_vec_s = np.linalg.eigh(L_s.toarray())

    L_st= L_s.transpose
    e_val_st,e_vec_st = np.linalg.eigh(L_st.toarray())
    
    return None #modify as needed



def modelA(G,x=0,i0=0.1,beta=1.0,gamma=1.0,tf=5,Nt=1000):
    """
    Question 2.2
    Simulate model A

    Input:
    G: Networkx graph
    x: node which is initially infected with i_x=i0
    i0: magnitude of initial condition
    beta,gamma: model parameters
    tf,Nt: Solutions are computed at Nt time steps from t=0 to t=tf (see code below)

    Output:
    iarray: N x Nt+1 Array containing i across network nodes at
                each time step.
    """

    N = G.number_of_nodes()
    iarray = np.zeros((N,Nt+1))
    tarray = np.linspace(0,tf,Nt+1)
    iarray[x,0] = i0
    A = nx.adjacency_matrix(G)
    

    def RHS(y,t,beta,gamma):
        """Compute RHS of modelA at time t
        input: y should be a size N array
        output: dy, also a size N array corresponding to dy/dt

        Discussion: add discussion here
        """
        dy = -beta*y + gamma* A@y *(1-y)
        
        return dy

    iarray = odeint(RHS,iarray[:,0],tarray,args=(beta,gamma),rtol=1e-4)
    
    iarray = np.transpose(iarray)
    
    return iarray

def modelB(G,x=0,i0=0.1,alpha=-0.01,tf=5,Nt=1000):
    A = nx.adjacency_matrix(G)
    Q = sparse.diags(A.toarray().sum(axis=1),format="csr")
    L = Q - A #generate L
    L=L.toarray() #try with sprase matrix, maybe faster
    N = G.number_of_nodes()
    
    #x=node0#fix this
    
    iarray = np.zeros(2*N) #long vector each half represenets one of s or i
    iarray[x] = i0
    tarray = np.linspace(0,tf,Nt+1) #time

    
    def RHS(y,t):
        #I
        #s
        dy = np.concatenate((alpha* L@y[N:], y[:N]), axis = None)
        return dy
    
    return np.transpose(odeint(RHS,iarray,tarray,rtol=1e-4)[:,N:]), np.transpose(odeint(RHS,iarray,tarray,rtol=1e-4)[:,:N]) #s and i equations


def transport(input=(None)):
    """Analyze transport processes (model A, model B, linear diffusion)
    on Barabasi-Albert graphs.
    Modify input and output as needed.
    """
    BAG_100 = nx.barabasi_albert_graph(100,5,seed=0)
    A100 = nx.adjacency_matrix(BAG_100) #need to find other instance of finding minimum node and change that to this because this method is better
    Q100 = A100.toarray().sum(axis=1)
    node0a=np.argmax(Q100)
    
    
    resultsa = modelA(BAG_100,node0a,i0=0.1,beta=0.5,gamma=0.1,tf=100)
    
    plt.figure(0)
    plt.plot(resultsa[node0a],label='initial node')
    plt.plot(resultsa[node0a +1],label='non initial node')
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelA')
    plt.legend(loc="lower right")


    plt.figure(1)
    plt.plot(np.transpose(resultsa))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelA for all nodes')
    
    resultsa = modelA(BAG_100,node0a,i0=0.1,beta=1,gamma=0.1,tf=100)
    plt.figure(2)
    plt.plot(np.transpose(resultsa))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelA for all nodes with larger beta')
    
    resultsa = modelA(BAG_100,node0a,i0=0.1,beta=0.5,gamma=0.5,tf=100)
    plt.figure(3)
    plt.plot(np.transpose(resultsa))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelA for all nodes with larger gamma')
      
    resultsb,resultsb2 = modelB(BAG_100,node0a,i0=0.1,tf=100)
    plt.figure(4)
    plt.plot(resultsb[node0a])
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelB')
    
    plt.figure(5)    
    plt.plot(np.transpose((resultsb)))
    plt.xlabel('time')
    plt.ylabel('flux')
    plt.title('flux against time for modelB')
    
    plt.figure(6)    
    plt.plot(np.transpose((resultsb2)))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.title('concentration against time for modelB')
  
    return None #modify as needed


if __name__=='__main__':
    #add code here to call diffusion and generate figures equivalent
    #to those you are submitting
    BAG_0 = nx.barabasi_albert_graph(2000,4,seed=0)
    print('2000 BAG generated')
    deg_dict = dict(BAG_0.degree)
    deg_dict_sort2 = {k: v for k, v in sorted(deg_dict.items(), key=lambda item: item[1],reverse=True)}
    node0 = list(deg_dict_sort2)[0]
    print('minimum node found')
    
    
    """___________""""
    #NOTE: run the next 4 rwgraph_analyze1 separately, otherwise plots are merged
    rwgraph_analyze1(8000,1,deg_plot=False)
    rwgraph_analyze1(8000,3,deg_plot=False)
    
    rwgraph_analyze1(200,200,deg_plot=False)
    
    rwgraph_analyze1(8000,200)

    G=None #modify as needed
    
    print('generating next batch of graphs')
    transport()
