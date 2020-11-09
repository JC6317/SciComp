""" CID:01063446
    Template code for part 2, contains 3 functions:
    codonToAA: returns amino acid corresponding to input amino acid
    DNAtoAA: to be completed for part 2.1
    pairSearch: to be completed for part 2.2
"""

import numpy as np
import timeit
import time
import matplotlib.pyplot as plt

def codonToAA(codon):
	"""Return amino acid corresponding to input codon.
	Assumes valid codon has been provided as input
	"_" is returned for valid codons that do not
	correspond to amino acids.
	"""
	table = {
		'ATA':'i', 'ATC':'i', 'ATT':'i', 'ATG':'m',
		'ACA':'t', 'ACC':'t', 'ACG':'t', 'ACT':'t',
		'AAC':'n', 'AAT':'n', 'AAA':'k', 'AAG':'k',
		'AGC':'s', 'AGT':'s', 'AGA':'r', 'AGG':'r',
		'CTA':'l', 'CTC':'l', 'CTG':'l', 'CTT':'l',
		'CCA':'p', 'CCC':'p', 'CCG':'p', 'CCT':'p',
		'CAC':'h', 'CAT':'h', 'CAA':'q', 'CAG':'q',
		'CGA':'r', 'CGC':'r', 'CGG':'r', 'CGT':'r',
		'GTA':'v', 'GTC':'v', 'GTG':'v', 'GTT':'v',
		'GCA':'a', 'GCC':'a', 'GCG':'a', 'GCT':'a',
		'GAC':'d', 'GAT':'d', 'GAA':'e', 'GAG':'e',
		'GGA':'g', 'GGC':'g', 'GGG':'g', 'GGT':'g',
		'TCA':'s', 'TCC':'s', 'TCG':'s', 'TCT':'s',
		'TTC':'f', 'TTT':'f', 'TTA':'l', 'TTG':'l',
		'TAC':'y', 'TAT':'y', 'TAA':'_', 'TAG':'_',
		'TGC':'c', 'TGT':'c', 'TGA':'_', 'TGG':'w',
	}
	return table[codon]


def DNAtoAA(S):
    """Convert genetic sequence contained in input string, S,
    into string of amino acids corresponding to the distinct
    amino acids found in S and listed in the order that
    they appear in S
    """    
    D = {} #empty dictionary
    
    for i in range(len(S)//3):
        if codonToAA(S[i:i+3]) not in D:
            D[codonToAA(S[i:i+3])] = S[i:i+3]
            
    AA = "".join(list(D))
    
    return AA


def generate_DNA(length): #length is a third of final length output
    DNA = "".join(np.random.choice(['G','T','A','C'],length*3)) #generate from GTAC and create multiple of 3
    return DNA

def DNAtoAA_time(l_start,l_end,n):
    
    results = np.zeros((n,2)) #empty list for results
    
    for i in range(n):
        l = l_start +  i * ((l_end-l_start)//n) #length
        S = generate_DNA(l) #generate DNA sequence
        t1 = time.time() #start time
        AA = DNAtoAA(S) 
        t2 = time.time() #end time
        results[i,0]= l
        results[i,1]= t2-t1
    
    plt.plot(results[:,0],results[:,1])
    plt.xlabel("sequence length")
    plt.ylabel("time (seconds)")
    return results
        
"""uncomment below to plot graph"""
#DNAtoAA_time(1000,50000,1500)
        
        

def char2base4(S):
    #convert string of GTAC into base 4 integer
    c2b = {} #dictionary
    c2b['A'] = 0
    c2b['C'] = 1
    c2b['G'] = 2
    c2b['T'] = 3
    L=[]
    for s in S: #for each character in string
        L.append(c2b[s]) #convert character to base 4. #?have to convert to string to add
    #L = int(''.join(L)) #join list of string integers and convert back to integer
    return L

def heval(L, base, prime):
    #convert list L to 'base' to base10 number mod prime
    f=0
    for l in L[:-1]:
        f = base*(l+f)
    h = (f + (L[-1])) % prime
    return h


def pairSearch(L,pairs):
    """Find locations within adjacent strings (contained in input list,L)
    that match k-mer pairs found in input list pairs. Each element of pairs
    is a 2-element tuple containing k-mer strings
    """
    
    locations = [] #empty list to store results
    #loop through for each tuple in the list of pairs
    #for pair in pairs:
    for i_pair in range(len(pairs)):
        pair = pairs[i_pair]
        #select each element of the chosen tuple and convert both into hash
        p1 = pair[0] #pattern 1
        p2 = pair[1] #pattern 2
        k = len(p1) #setting k length
        
        Ip1 = char2base4(p1) #base 4 representation
        Ip2 = char2base4(p2)
        hp1 = heval(Ip1,4,97) #base 10 hashpattern1, prime 
        hp2 = heval(Ip2,4,97)
        
        #loop through each DNA sequence string in list L
        for i_dna in range(len(L)-1): #-1 because we do not need to check thebottom element
            dna = L[i_dna]
            #rolling hash on string with pair1
            #for roll in range(len(dna)-2): #rolling index for STARTING LOCATION of 3 letter string
            for roll in range(len(dna)-k+1): #rolling index for STARTING LOCATION of 3 letter string
                rollhash = heval(char2base4(dna[roll:roll+k]),4,97) #rolling hash of 3 letter stringin dna sequence
                if hp1 == rollhash: #if hash pattern 1 matches rolling hash
                    if p1 == dna[roll:roll+k]: #check for hash collision
                        belowstringhash = heval(char2base4(L[i_dna+1][roll:roll+k]),4,97) #calculate hash for below string (L[i+1])
                        if hp2 == belowstringhash:
                            if p2==L[i_dna+1][roll:roll+k]:
                                #save
                                match = [roll,i_dna,i_pair] #create 3 element sublist 
                                locations.append(match) #append


    return locations

"""Testing pairsearch"""
L = ["GCAATTCGT","TCGTTGATC", "GTCAAGTAC"]
pairs = [("TCG","GAT"),("ATC","TAC")]
pairSearch(L,pairs)