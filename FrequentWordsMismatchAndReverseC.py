
import math
import numpy as np

#Frequent words with mismatches and reverse complements, takes in dna strand, length of k-mer and # of mismatches allowed
def FrequentWordsWithMismatchesAndReverseComplements(dna, k, d) : 
    dna = dna.replace("\n","")
    #define the array with each equal to 0
    freq_map = {}
    
    for i in range((4**k)) :
        freq_map[i] = 0
    #now loop through the DNA strand and map them
    for sub in [dna[i: i + k] for i in range(len(dna) -(k -1))] :
        #for each sub get the neighbors of it with d mismatches
        close = neighbors(sub, d)
        for c in close : 
            #get the numValue of c
            num = patternToNumber(c)
            freq_map[num] = freq_map.get(num) + 1   
    return mostFrequentWithCompliment(freq_map, k)

#This will go through our freq array and get the total occurences of pattern + total occurences of compliment(pattern)
def mostFrequentWithCompliment(freq_arr, k) :
    
    max = 0
    most_freq = []
    for key in freq_arr :
        pat = NumberToPattern(key,k)
        #get the reverse compliment
        rev = compliment(pat)
        revIndex = patternToNumber(rev)
        
        if (freq_arr.get(key) + freq_arr.get(revIndex)) > max :
            max = freq_arr.get(key) + freq_arr.get(revIndex)
            most_freq = []
            most_freq.append(NumberToPattern(key,k))
        elif (freq_arr.get(key) + freq_arr.get(revIndex)) == max :
            most_freq.append(NumberToPattern(key,k))
            
    return most_freq

#converts the dna from its lexographic index to pattern
#takes in the index and the number of bases in the pattern
def NumberToPattern(index, num_bases) :
    
    order = {0 : 'A' , 1 : 'C' , 2 : 'G' , 3 : 'T'}
    
    pattern = ""
    
    for i in range(num_bases) : 
        
        power = 4 ** (num_bases-i-1)
        
        temp = math.floor(index / power)

        pattern += str(order.get(temp))
        
        index = index - (power * temp)
    return pattern

#calculates the lexographic index of the DNA pattern
def patternToNumber (pattern) :   
    index = 0    
    order = {'A' : 0 , 'C' : 1 , 'G' : 2 , 'T' : 3}    
    i = len(pattern) #the position we are in the pattern
    
    #scan through the pattern  
    for base in pattern :
        index = index + int(order.get(base)) * (4**(i-1))
        i = i -1
    
    return index

#returns a list of the mismatches assocaiated with a pattern with up to d mismatches
def neighbors(pattern, d):
    close = [pattern]
    for i in range(d) :
        #loop through what is in close
        res = []
        for pat in close :
            for neigh in changeOne(pat) :
                res.append(neigh)
        for r in res :
            close.append(r)
    close = np.unique(close)
    return close


#returns all the mismatches when just one change occurs 
def changeOne(sub): 
    replace = {'A' :['T','C','G'],'T' :['A','C','G'], 'G' :['A','C','T'] ,'C' :['A','T','G'] }
    res = []
    index = 0
    pattern = list(sub)
    for i in pattern :
        #replace the nucleotide at index with the three other options
        for j in replace.get(i): 
            temp = list(sub) 
            temp[index] = j 
            res.append("".join(temp))
        index +=1
    return(res)


#returns the compliment of the DNA strand
def compliment(dna): 
    comp = ''
    map = {"A":"T","T":"A","G":"C","C":"G"}
    for c in reversed(dna):
        comp += str(map.get(c))
        
    return comp


