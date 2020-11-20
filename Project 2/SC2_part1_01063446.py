"""Scientific Computation Project 2, part 1
Your CID here: 01063446
"""
import numpy as np
from collections import deque

def flightLegs(Alist,start,dest):
    """
    Question 1.1
    Find the minimum number of flights required to travel between start and dest,
    and  determine the number of distinct routes between start and dest which
    require this minimum number of flights.
    Input:
        Alist: Adjacency list for airport network
        start, dest: starting and destination airports for journey.
        Airports are numbered from 0 to N-1 (inclusive) where N = len(Alist)
    Output:
        Flights: 2-element list containing:
        [the min number of flights between start and dest, the number of distinct
        jouneys with this min number]
        Return an empty list if no journey exist which connect start and dest
    """
    Alen = len(Alist) #length of Alist
    index = list(range(Alen)) #list of airport indexes
    Adict = dict(zip(index,Alist)) #zip into dictionary
    #Adict = dict.fromkeys(range(Alen),large value) possibly qiucker?
    
    Flights = [] #output

    checked = [] #empty list
    tocheck = deque() #empty deque for a queue of paths to be checked
    tocheck.append([start])
    
    if start == dest: #trivial case
        return "trivial case explored"
    
    while tocheck: #while there are elements in tocheck
        path = tocheck.popleft() #pop out a path to be checked
        node = path[-1] #extract last element of path, this node represents an airport
        if node not in checked: #check if this path terminates in a new node/airport
            adjacents = Adict[node] #look up in dictionary airport index and return list of adjacent airports
            for n in adjacents: #extracts nodes inside list of adjacents
                if n not in checked:
                    new_path = path +[n] #add neighbour to this current path
                    tocheck.append(new_path) #put this new path into the queue to be checked
                    if n == dest: #if this neighbour is the final destination, then the new path is the full path
                        if len(Flights)==0:
                            Flights.append(new_path)
                        elif len(new_path)==len(Flights[-1]):
                            Flights.append(new_path)
                        else:
                            return (len(Flights[0]),len(Flights))
                        
            checked.append(node)

    return Flights #finds all paths need to modify to stop after a path of non shortest length is found


def safeJourney(Alist,s,d):
    """
    Question 1.2 i)
    Find safest journey from station start to dest
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the density for the connection.
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing safest journey and safety factor for safest journey
    """
    #Initialize dictionaries
    dinit = 10**6
    Edict = {} #Explored nodes
    Udict = {} #Unexplored nodes
    path = [[] for l in Alist]

    Alen = len(Alist) #length of Alist
    dinits = [dinit]*Alen #list of airport indexes
    Udict = dict(zip(list(range(Alen)),dinits)) #zip into dictionary
    Udict[s] = 0
    path[s] = [s]
        
    #Main search
    while len(Udict)>0:
        #Find node with min d in Udict and move to Edict
        dmin = dinit
        for n,w in Udict.items():
            if w<dmin:
                dmin=w
                nmin=n
        Edict[nmin] = Udict.pop(nmin)
        print("moved node", nmin)

        #Update provisional distances for unexplored neighbors of nmin
        
        #for n,w in G.adj[nmin].items():
        for item in Alist[nmin]: #nminth element is a list of two element tuples (node, weight)
            n = item[0] #first elt of tuple is node/neighbour
            w = item[1] #2nd elt is density/weigh
            #for n,w in etc_______________________-
            
            if n in Edict:
                pass
            elif n in Udict:
                #dcomp = dmin + w
                dcomp = max(w,dmin) #take largest value to record most dangerous segment
                if dcomp<Udict[n]:
                    print(Udict)
                    Udict[n]=dcomp
                    path[n] = path[nmin] + [n]
                    #path[n].extend(path[nmin])
                    #path[n] = path[nmin]
                    
                    #path[n].append(n) #n not nmin
                    print(path)
           # else:
                #dcomp = dmin + w
            #    dcomp = max(w,dmin)
             #   Udict[n] = dcomp
                #path[n].extend(path[nmin])
                #path[n].append(nmin) 
               
            if nmin == d: #if current node is destination
                return path[d],Edict[d]
    return [] #no path



def shortJourney(Alist,s,d):
    """
    Question 1.2 ii)
    Find shortest journey from station start to dest. If multiple shortest journeys
    exist, select journey which goes through the smallest number of stations.
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the time for the connection (rounded to the nearest minute).
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing shortest journey and duration of shortest journey
    """
    """Find shortest distances to s in weighted graph, G"""
        
    #Initialize dictionaries
    dinit = 10**6
    Edict = {} #Explored nodes
    Udict = {} #Unexplored nodes
    path = [[] for l in Alist]

    Alen = len(Alist) #length of Alist
    dinits = [dinit]*Alen #list of airport indexes
    Udict = dict(zip(list(range(Alen)),dinits)) #zip into dictionary
    Udict[s] = 0
    path[s] = [s]
    
    #Main search
    while len(Udict)>0:
        #Find node with min d in Udict and move to Edict
        dmin = dinit
        for n,w in Udict.items():
            if w<dmin:
                dmin=w
                nmin=n
        Edict[nmin] = Udict.pop(nmin)
        print("moved node", nmin)

        #Update provisional distances for unexplored neighbors of nmin     
        for item in Alist[nmin]: #nminth element is a list of two element tuples (node, weight)
            n = item[0] #first elt of tuple is node/neighbour
            w = item[1] #2nd elt is density/weigh
            #for n,w in etc_______________________-
            
            if n in Edict:
                pass
            elif n in Udict:
                #key difference below
                dcomp = (w+dmin) #take sum as you go along
                if dcomp<Udict[n]:
                    print(Udict)
                    Udict[n]=dcomp
                    path[n] = path[nmin] + [n]
                    print(path)           
            if nmin == d: #if current node is destination
                return [path[d],Edict[d]]
    return [] #no path
    
def cheapCycling(SList,CList):
    """
    Question 1.3
    Find first and last stations for cheapest cycling trip
    Input:
        Slist: list whose ith element contains cheapest fare for arrival at and
        return from the ith station (stored in a 2-element list or tuple)
        Clist: list whose ith element contains a list of stations which can be
        cycled to directly from station i
    Stations are numbered from 0 to N-1 with N = len(Slist) = len(Clist)
    Output:
        stations: two-element list containing first and last stations of journey
    """
    N = len(CList)
    nodes= list(range(N))
    Udict = dict(zip(nodes,CList)) #dictionary nodes:neighbours
    queue = deque() #things to check
    checked = {} #checked nodes to avoid double checking
    stations=np.zeros((N,2))#node:arrival output

    min_arr = 1000000 #initialise fat value
    arr_node=-1
    
    min_dep = 1000000
    dep_node=-1

    while Udict:
        #node = Udict.pop(next(iter(Udict))) #extract first node
        node = next(iter(Udict)) #select first node
        queue.append(node) #add node to queue
        while queue: #while queue is non empty
            node = queue.popleft() #set/extract node to element of queue
            Udict.pop(node) #make sure is also removed from overarching dict
            for nb in CList[node]: #neighbours of node
                if nb not in checked:
                    if SList[nb][0] <min_arr: #check if new minimum
                        min_arr = SList[nb][0]
                        arr_node = nb
                    if SList[nb][1] <min_dep: #same but for departure
                        min_dep = SList[nb][1]
                        dep_node= nb
                    queue.append(nb)
                    checked[nb] = 1
            checked[node]=1
                    
        stations[list(checked.keys()),0] = arr_node #dropping in the cheapest arr and dep nodes for all nodes in connected part
        stations[list(checked.keys()),1] = dep_node
        checked={} #reset checked to empty for new connected part IMPORTANT!
        min_arr=1000000 #reset minimum values
        min_dep=1000000 #reset min dep values

    return stations





if __name__=='__main__':
    #add code here if/as desired
    L=None #modify as needed
