# == This module of functions mainly used for z-score calculation
import random
import itertools
import numpy as np
import buildMatrix as bm
import utility
import time

AA_LIST = list('ARNDCEQGHILKMFPSTWYV')

class CheckZScore(object):
    
    # initialize the matrix in the data structure of dict of list
    def __init__(self, MAX_LENGTH, FGset, fixed_pairs, NT_BG_freq_matrix, CT_BG_freq_matrix, original_log_bp,
                 NT_pseudo_count_matrix, CT_pseudo_count_matrix, log_factorial_map):
        self.MAX_LENGTH = MAX_LENGTH
        self.fixed_pairs = fixed_pairs
        self.NT_BG_freq_matrix = NT_BG_freq_matrix
        self.CT_BG_freq_matrix = CT_BG_freq_matrix
        self.original_log_bp = original_log_bp
        self.NT_pseudo_count_matrix = NT_pseudo_count_matrix
        self.CT_pseudo_count_matrix = CT_pseudo_count_matrix
        self.log_factorial_map = log_factorial_map 
        self.FGset = FGset # alias_fg_set
        self.FG_list_char = list() # list of char list of FGset
        for pep in FGset:
            self.FG_list_char.append(list(pep))
        #Initialize the list of groups of different length for faster shuffle
        self.groups_FG_list_char = {}
        self.key_length = []    
        #Group peptides into groups based on their length
        sorted_FG_list_char = sorted(self.FG_list_char, key=lambda a : len(a))
        for k, g in itertools.groupby(sorted_FG_list_char, lambda a:len(a)):
            group = np.array(list(g))
            self.groups_FG_list_char[k] = group      # Store group iterator as a list
            self.key_length.append(k)
        
        
            
            
    # == Function Summary
    #   Calculate Z score based on randomization of FG 
    #   then calculate mean and standard deviation for the randomized FG
    #   Parameters: 
    #   -- pep_fg, peptides in foreground
    #   -- residue, amino acid in consideration
    #   __ position, 
    #   -- orig_count, 
    #   Return:
    #   function return a shuffled matrix
    #   example: the shuffled matrix looks like 
    
     
    def calculate_ZScore(self):
        
        #count iteration
        itr = 0
        bp_list = []
        start_time = time.time()
        mean = 0 # average of bp_list so far
        most_sig_bp = -100 # initialize a significant number for converge calculation
        
        while (True):
            if (itr > 49):
                break
            self.shuffle_matrix() # shuffle matrix for calculation
            fg_set_shuffled = self.convert_to_dataset(self.groups_FG_list_char)
            #print "fg_set_shuffled=", fg_set_shuffled
            most_sig_bp = self.searchMostSignificantBP(fg_set_shuffled)
            #print("sig_bp = "), most_sig_bp
            
            bp_list.append(most_sig_bp)
            # When a random bp less than the one to compare, show there is no significance
#            if most_sig_bp < self.original_log_bp:
#                return 0
            if itr > 1:
                mean = sum(bp_list)/len(bp_list)
            if (itr == 4):
                differences = [(value - mean)**2 for value in bp_list]
                mean_differences = sum(differences)/len(differences)
                stdev = (mean_differences) ** 0.5
                #print("stdev = "), stdev
                # Calculate z-score
                if stdev != 0:
                    z_score = abs(self.original_log_bp - mean)/stdev
                    if (z_score > 10):
                        break
            itr += 1
            #print("itr = "), itr
            
            
        # Calculate mean and standard divation
        end_time = time.time()
        runtime = (end_time - start_time)/60
        #print("one cycle of shuffle = "), runtime, " min"
        
        differences = [(value - mean)**2 for value in bp_list]
        #print "differences = ", differences
        mean_differences = sum(differences)/len(differences)
        stdev = (mean_differences) ** 0.5
        print ("stdev = "), stdev
        print ("bp_list = "), bp_list
        # Calculate z-score
        if stdev != 0:
            z_score = abs(self.original_log_bp - mean)/stdev
        else:
            z_score = 1000000
        
        return z_score
    
    
     # == Function Summary
    #   Shuffle matrix, anchor position can be designated and kept unchanged 
    #   used for randomization and permutation
    #   Parameters: 
    #   -- pos, an list storing all positions will be unchanged in data structure set
    #   -- shuffled_matrix, peptides saved in the list of list data structure 
    #   Return:
    #   function return a shuffled matrix
    #   example: the shuffled matrix looks like 
    #   array([[ 'A', 'C', 'D', 'E', 'S', 'R', 'L'],
    #          [ 'F', 'R', 'D', 'E', 'S', 'R', 'L', 'Y'],
    #          [ 'T', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W'],
    #          [ 'S', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W'],
    #          [ 'V', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W']]) 
     
    def shuffle_matrix(self):
        
        #Put all positions in fixed_pairs into a set
        fixed_pos_set = set([])
        for pair in self.fixed_pairs:            
            fixed_pos_set.add(pair[1])       
        
        #print "fixed_pos_set = ", fixed_pos_set
        # Before shuffle, construct a long list to store all shuffled aa
        shuffle_peps = []
        # debug printing
        #print ("groups_FG_list_char = "), self.groups_FG_list_char
        
        # Use numpy to construct the array for shuffle
        for key in self.key_length:
            for p in range(key):
                if (p not in fixed_pos_set) and ((p-key) not in fixed_pos_set):
                    shuffle_peps.extend(self.groups_FG_list_char[key][:,p])
                
        # Shuffle this list in place (without replacement)
        random.shuffle(shuffle_peps)  
        #print "shuffle_peps = ", len(shuffle_peps)
        
        # Use numpy to reconstruct data after shuffle
        start = 0 #The iteration numbers that processed
        for key in self.key_length:
            rows = len(self.groups_FG_list_char[key])
            for p in range(key):
                if (p not in fixed_pos_set) and ((p-key) not in fixed_pos_set):
                    self.groups_FG_list_char[key][:,p] = shuffle_peps[start:start+rows]
                    start += rows

    
    # == Function Summary
    #   Convert foreground and background data set into matrix, 
    #   used for randomization and permutation
    #   Parameters: 
    #   -- dataset, peptides saved in the set([]) data structure 
    #   Return:
    #   function return a matrix, i.e. list of list with one letter as an element
    #   example: the matrix looks like 
    #   list( [[ 'A', 'C', 'D', 'E', 'S', 'R', 'L'],
    #          [ 'F', 'R', 'D', 'E', 'S', 'R', 'L', 'Y'],
    #          [ 'T', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W'],
    #          [ 'S', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W', 'U'],
    #          [ 'V', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W', 'U', 'V']]) 
    
    def charMatrixAlignedNT(self, seqset):
        pepArray_str = np.array(list(seqset), dtype=str)
        pepArray_NT_str = np.core.defchararray.ljust(pepArray_str, self.MAX_LENGTH)
        pepArray_NT_char = pepArray_NT_str.view('S1').reshape(pepArray_NT_str.size, -1)
        return pepArray_NT_char

    # generate backward aligned matrix
    def charMatrixAlignedCT(self, seqset):
        pepArray_str = np.array(list(seqset), dtype=str)
        pepArray_CT_str = np.core.defchararray.rjust(pepArray_str, self.MAX_LENGTH)
        pepArray_CT_char = pepArray_CT_str.view('S1').reshape(pepArray_CT_str.size, -1)
        return pepArray_CT_char
      
   
    
    # == Function Summary
    #   Convert foreground and background data matrix into set datastructure, 
    #   used for compatible to other utility methods
    #   Parameters: 
    #   -- matrix saved in the list of list data structure 
    #   example: the matrix looks like 
    #   list( [[ 'A', 'C', 'D', 'E', 'S', 'R', 'L', 'I'],
    #          [ 'F', 'R', 'D', 'E', 'S', 'R', 'L', 'Y', 'W'],
    #          [ 'T', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W', 'U'],
    #          [ 'S', 'C', 'D', 'E', 'S', 'R', 'L', 'I', 'W', 'U', 'V'],
    #          [ 'V', 'C', 'D', 'E', 'S', 'R', 'L']]) 
    #   Return:
    #   function return a set structure
    
     
    def convert_to_dataset(self, data_list_char):
        pep_dataset = set([])
        for key in data_list_char:
            for pep_list in data_list_char[key]:
                pep_string = "".join(pep_list)
            #print("pep_list = "), pep_list
                pep_dataset.add(pep_string)   
            #print("pep_matrix = "), pep_matrix
        return pep_dataset
                
    # == Function Summary
    #   Calculate the time observe r/p pair in a matrix
    #   Parameters: 
    #   -- fg, shuffled foreground, a list of list
    #   -- bg, background without shuffling
    #   Return:
    #   function return the most significant binomial probability with the current foreground
    
     
    def searchMostSignificantBP(self, shuffled_set):
        most_sig_bp = 0
        
        #Put all positions in fixed_pairs into a set
        fixed_positions = set([])
        for pair in self.fixed_pairs:
            fixed_positions.add(pair[1])
        
        matrix = bm.Matrix(self.MAX_LENGTH, shuffled_set)
        NT_shuffled_charMatrix = self.charMatrixAlignedNT(shuffled_set)
        CT_shuffled_charMatrix = self.charMatrixAlignedCT(shuffled_set)
        NT_shuff_count_matrix = matrix.aaCountMatrix(NT_shuffled_charMatrix)
        CT_shuff_count_matrix = matrix.aaCountMatrix(CT_shuffled_charMatrix)  
         
        
        # search the most significant binomial probability
        for aa in AA_LIST:
            for pos in range(self.MAX_LENGTH):
                if pos not in fixed_positions:
                    k = NT_shuff_count_matrix[aa][pos][0]
                    n = NT_shuff_count_matrix[aa][pos][1]
                    p = self.NT_BG_freq_matrix[aa][pos]
                    # if foreground completely deplete one rp pair, set the binomial probability to 0 (its log form)
                    if k == 0 or p == 0:
                        log_bp = 0
                    elif k == n and p == 1:
                        log_bp = 0    
                    else:
                        log_bp = utility.binomial_log_prob(n, k, p, self.log_factorial_map)
                        
                    #print "this_log_bp =", log_bp
                    
                    if log_bp < most_sig_bp:
                        most_sig_bp = log_bp
                        
        # search the most significant binomial probability
        for aa in AA_LIST:
            for pos in range(self.MAX_LENGTH):
                if pos not in fixed_positions:
                    k = CT_shuff_count_matrix[aa][pos][0]
                    n = CT_shuff_count_matrix[aa][pos][1]
                    p = self.CT_BG_freq_matrix[aa][pos]
                    # if foreground completely deplete one rp pair, set the binomial probability to 0 (its log form)
                    if k == 0 or p == 0:
                        log_bp = 0
                    elif k == n and p == 1:
                        log_bp = 0    
                    else:
                        log_bp = utility.binomial_log_prob(n, k, p, self.log_factorial_map)
                        
                    #print "this_log_bp =", log_bp
                    if log_bp < most_sig_bp:
                        most_sig_bp = log_bp
        
        return most_sig_bp
    
                

            
