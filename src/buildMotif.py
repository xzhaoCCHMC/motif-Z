# == Summary
#   This module perform iterative motif finding algorithm and produce complete motif list for a input data set 
#   motifs were compiled in regular expression format, order of pairs saved in a list and passed to webpage to process
#
# 


import generateOneMotif as gom
import utility
from math import exp
import copy

# Global variables
global debug
debug = False

# Constant and default values
PSEUDO_N_SIZE = 50
AA_LIST = list('ARNDCEQGHILKMFPSTWYV')


class MotifSet:
    
    #   -- MAX_LENGTH, the length of peptide
    #   -- BINOMIAL_PROB_THRESHOLD, value used as threshold to identify significant pair
    def __init__(self, max_length, binomial_prob_threshold, NT_pseudo_count_matrix, CT_pseudo_count_matrix, 
                 threshold_fg, original_fg, original_bg, min_length, log_factorial_map, medium_pep_length, ZSCORE_FLAG):
        self.motifList = list() # rp pairs in list
        self.motifSetRef = {} # Dictionary to store result and parameter information
        self.motifDict = {} # match motif with pairs
        self.motifid = 0
        self.rCount = 0
        self.Flag = 0
        self.MAX_LENGTH = max_length
        self.MINIMUM_LENGTH = min_length
        self.ZSCORE_FLAG = ZSCORE_FLAG
        self.BINOMIAL_PROB_THRESHOLD = binomial_prob_threshold
        self.NT_pseudo_count_matrix = NT_pseudo_count_matrix # pseudo count matrix from N-terminus alignment
        self.CT_pseudo_count_matrix = CT_pseudo_count_matrix # pseudo count matrix from C-terminus alignment
        self.threshold_fg = threshold_fg  # User defined minimum peptide sequences in foreground to make the cut
        self.original_fg = original_fg # Original foreground in data structure set
        self.original_bg = original_bg # Original background in data structure set
        self.log_factorial_map = log_factorial_map # Precompute factorial values 
        self.medium_pep_length = medium_pep_length # The majority of peptide length in dataset
	self.rpSeqList = list() # list for order of sequential rp in the motifs for presentation
                
    # == Function Summary
    #   Perform iterative motif finding algorithm on input dataset and output motif set derived
    #   recursion will update both foreground and background based on each motif identified in previous search circle 
    
    def produce_all_motif(self, foreground_peps_set, background_peps_set):
        
        if self.rCount == 0:
            print("####################################################################################")
            print("Round " + str(self.rCount) + ": "),
            print("FG" + str(self.rCount) +" = "+ str(len(foreground_peps_set)) + "\t BG" + str(self.rCount) +
               " = " + str(len(background_peps_set)))
        
#        # for debug purpose, stop at 10 recursion steps
#        if self.rCount == 30:
#            print("for debugging, stop here")
#            return self.motifSet, self.motifSetRef
        
        
        # Loop over until no peptide motif can be identified  
        """if self.Flag == 1:
            print("Note: All significant residue and position pairs have been found for the input data file")
            return self.motifList, self.motifSetRef"""
        if len(background_peps_set) == 0:
            print("Warning: Empty background, suspend search!")
            return self.motifList, self.motifSetRef, self.motifDict
        elif len(foreground_peps_set) == 0:
            print("Warning: Empty foreground, suspend search!")
            return self.motifList, self.motifSetRef, self.motifDict
        
        # This variables just for debugging
        start_fg = foreground_peps_set
        start_bg = background_peps_set
        
        # Initialize two alias set for replace fixed rp with 'X'
        alias_fg_set = copy.deepcopy(foreground_peps_set)
        alias_bg_set = copy.deepcopy(background_peps_set)
        
        # initialize fg_branch and bg_branch
        fg_branch = foreground_peps_set
        bg_branch = background_peps_set
                
        
        # Produce one motif at a time
        # This used for debugging
        n_rp = 0 # counting the number of residue/position pair found
        # parameter for output of binomial probability for the last residue/position pair, initialize to 1
        last_p_sig_pair = 1
        # initialize the score for each motif
        score = 0.0
        l_fold = 0
        
        
        motif = gom.Motif(self.MAX_LENGTH, self.BINOMIAL_PROB_THRESHOLD, self.NT_pseudo_count_matrix, 
                          self.NT_pseudo_count_matrix, self.threshold_fg, self.original_fg, self.original_bg, 
                          self.MINIMUM_LENGTH, self.log_factorial_map, self.medium_pep_length, self.ZSCORE_FLAG) 
                       
        [motif_identified, searched_fg_set, searched_bg_set, debug_n, fg_vs_branch, bg_vs_branch, binomial_p, 
         l_enrichment, score, Stop] = motif.produce_one_motif(fg_branch, foreground_peps_set, bg_branch, 
                background_peps_set, self.Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set)  
         
        # Check if all motifs have been exhausted, if so, stop recursion    
        if Stop == 1 or not motif_identified: 
            print("Note: All significant pairs identified and search is exhaustive for the input data file")  
            return self.motifList, self.motifSetRef, self.motifDict
        
        if motif_identified in self.motifList:
            print("Note: motif start to repeat!")  
            return self.motifList, self.motifSetRef, self.motifDict
        
        # string for final regular expression
        rexp_motif_str = utility.format_to_rexp(motif_identified, self.MAX_LENGTH, self.MINIMUM_LENGTH)
 	
	# list of each residue position in the rexp and the corresponding iteration number
	rpSeqList = utility.rp_iteration_id(motif_identified, rexp_motif_str)
        #regular expression of the motif from String format
        '''motif_rexp_format = re.compile(rexp_motif_str)
        fold_enrichment, num_in_fg, num_in_bg = utility.fold_enriched(self.original_fg, 
                                                                      self.original_bg, motif_rexp_format, "rexp")'''
 
        fold_enrichment, num_in_fg, num_in_bg = utility.fold_enriched(self.original_fg, 
                                                                      self.original_bg, motif_identified, "pairs")
        
        #Calculate p value for the motif, using binomial probability for the whole motif
        n = len(self.original_fg)
        k = num_in_fg
        p = float(num_in_bg)/len(self.original_bg)
        log_bp_motif = utility.binomial_log_prob(n, k, p, self.log_factorial_map)
        
        self.motifSetRef[rexp_motif_str] = [fg_vs_branch, bg_vs_branch, binomial_p, l_enrichment, score, 
                                        fold_enrichment, num_in_fg, num_in_bg, log_bp_motif]
        self.motifDict[rexp_motif_str] = [motif_identified, rpSeqList]
        
        
        # print out information for the debug file to track the fg and bg size
        print("Final Motif = "), rexp_motif_str
        print("rp pairs in motif = "), motif_identified
	print("rpSeqList="), rpSeqList
        print("(" + str(len(searched_fg_set)) + "/ " + str(len(searched_bg_set))) + ") \t",
        motif_enrichment, cur_fg_size, cur_bg_size = utility.fold_enriched(start_fg, start_bg, motif_identified, "pairs")
        print("Fold enrichment in this round = "),
        if motif_enrichment == 'NA':
            print('NA')
        else:
            print("{0:.2f}".format(motif_enrichment))
        
            
        #update the motif
        if motif_identified:
            self.motifList.append(rexp_motif_str)
            
        
        print("Total Motif = "),
        print(self.motifList),
        
        print("\t Total Motif size = "),
        print(len(self.motifList))
        
        #Save the current foreground and background set as fg_branch and bg_branch
        fg_branch = foreground_peps_set
        bg_branch = background_peps_set
        
        #update the foreground and background sets by deducting the found motif patterns
        foreground_peps_set = foreground_peps_set - searched_fg_set
        background_peps_set = background_peps_set - searched_bg_set
        
        self.rCount += 1 # Recursion counting
        self.motifid += 1 # use motifid as key for dictionary of motifSetRef
        print("#################################################################################################")
        print("Round " + str(self.rCount) + ": "),
        print("FG" + str(self.rCount) + " = FG" + str(self.rCount - 1) + " - " + " FG" + str(debug_n) + " = " 
              + str(len(foreground_peps_set)) + "\t BG" + str(self.rCount) + " = " + str(len(background_peps_set)))
        
    
        return MotifSet.produce_all_motif(self, foreground_peps_set, background_peps_set)    
        
        
        
    
   
        
