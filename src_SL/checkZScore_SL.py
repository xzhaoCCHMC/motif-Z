# == This module of functions mainly used for z-score calculation
import random
import itertools
import numpy as np
import buildMatrix_SL as bm
import utility
import time
#from collections import Counter
from python27_modules import Counter

AA_LIST = list('ARNDCEQGHILKMFPSTWYV')

class CheckZScore(object):
    
    # initialize the matrix in the data structure of dict of list
    def __init__(self, PEP_LENGTH, FGset, fixed_pairs, BG_freq_matrix, original_log_bp,
                 log_factorial_map):
        self.PEP_LENGTH = PEP_LENGTH
        self.fixed_pairs = fixed_pairs
        self.BG_freq_matrix = BG_freq_matrix
        self.original_log_bp = original_log_bp
        self.log_factorial_map = log_factorial_map 
        self.FGset = FGset # alias_fg_set
        self.total_peptides=len(FGset) # number of peptides in the FG
        self.FG_list_char = list() # list of char list of FGset
        for pep in FGset:
            self.FG_list_char.append(list(pep))
        self.FG_array_char = np.array(self.FG_list_char)
            
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
        most_sig_bp = -100 # initialize a significant number for converge calculation
        
        while (True):
            if (itr > 49):
                break
            
            self.shuffle_matrix() # shuffle matrix for calculation
            #fg_set_shuffled = self.convert_to_dataset(self.FG_list_char)
            #print "fg_set_shuffled=", self.FG_array_char
            #most_sig_bp = self.searchMostSignificantBP(fg_set_shuffled)
            most_sig_bp = self.searchMostSignificantBP()
            #print("sig_bp = "), most_sig_bp
            
            bp_list.append(most_sig_bp)
            
#            if itr == 0:
#                mean_prev = most_sig_bp# initialize mean

            #if most_sig_bp < self.original_log_bp:
                #return 0
            
#                
            if itr > 1:
                mean=np.mean(bp_list)
                print "at round " + str(itr) + ", mean = ", mean
#                print "previous mean = ", mean_prev
#                if abs(mean-mean_prev) < 0.01:
#                    break;
#                else:
#                    mean_prev = mean
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
    #          [ 'F', 'R', 'D', 'E', 'S', 'R', 'L'],
    #          [ 'T', 'C', 'D', 'E', 'S', 'R', 'L'],
    #          [ 'S', 'C', 'D', 'E', 'S', 'R', 'L'],
    #          [ 'V', 'C', 'D', 'E', 'S', 'R', 'L']]) 
     
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
        for p in range(self.PEP_LENGTH):
            if p not in fixed_pos_set:
                shuffle_peps.extend(self.FG_array_char[:,p])
                
        # Shuffle this list in place (without replacement)
        random.shuffle(shuffle_peps)  
        #print "shuffle_peps = ", len(shuffle_peps)
        
        # Use numpy to reconstruct data after shuffle
        start = 0 #The iteration numbers that processed
        for p in range(self.PEP_LENGTH):
            if p not in fixed_pos_set:
                self.FG_array_char[:,p] = shuffle_peps[start:start+self.total_peptides]
                start += self.total_peptides
        
        
    
    # == Function Summary
    #   Convert foreground and background data matrix into set datastructure, 
    #   used for compatible to other utility methods
    #   Parameters: 
    #   -- matrix saved in the list of list data structure 
    #   example: the matrix looks like 
    #   list( [[ 'A', 'C', 'D', 'E', 'S', 'R', 'L', 'I'],
    #          [ 'F', 'R', 'D', 'E', 'S', 'R', 'L', 'Y'],
    #          [ 'T', 'C', 'D', 'E', 'S', 'R', 'L', 'I'],
    #          [ 'S', 'C', 'D', 'E', 'S', 'R', 'L', 'I'],
    #          [ 'V', 'C', 'D', 'E', 'S', 'R', 'L', 'Y']]) 
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
    
     
    def searchMostSignificantBP(self):
        most_sig_bp = 0
        
        #Put all positions in fixed_pairs into a set
        fixed_positions = set([])
        for pair in self.fixed_pairs:
            fixed_positions.add(pair[1])
        
        count_matrix_dict = dict()
        for aa in AA_LIST:
            count_matrix_dict[aa] = list()
            
        for col in range(self.PEP_LENGTH):
            count_curr_col = Counter.Counter(self.FG_array_char[:, col])
            if count_curr_col.has_key(' '):
                del count_curr_col[' ']
            if count_curr_col.has_key('X'):
                del count_curr_col['X']
            
            for aa in AA_LIST:
                if count_curr_col.has_key(aa):
                    count_aa = count_curr_col[aa]
                else:
                    count_aa = 0
                count_matrix_dict[aa].append(count_aa)           
            
        shuff_count_matrix = count_matrix_dict
        
         
        # search the most significant binomial probability
        for aa in AA_LIST:
            for pos in range(self.PEP_LENGTH):
                if pos not in fixed_positions:
                    k = shuff_count_matrix[aa][pos]
                    n = self.total_peptides
                    p = self.BG_freq_matrix[aa][pos]
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
    
