'''
Created on Mar 24, 2013
# == Summary
#   This module contains the major iteration and recursion procedure for Motif-Z algorithm
#   Perform iterative motif finding algorithm on input file and writes output
#   This function module returns a motif as rp pairs in a list, with rp in a sequential positions
#   
@author: zhaoxh
'''
        
import findSignificantPair as fsp  
import utility
import re
from math import exp
import time

# Global variables
global debug
debug = False

# Constant and default values
PSEUDO_N_SIZE = 50
AA_LIST = list('ARNDCEQGHILKMFPSTWYV')


class Motif(object):
    '''
    = Summary
    The class is to build a single motif based on significant pairs identified
    -- MAX_LENGTH, the length of peptide
    -- BINOMIAL_PROB_THRESHOLD, value used as threshold to identify significant pair
    '''
    
    def __init__(self, max_length, binomial_prob_threshold, NT_pseudo_count_matrix, CT_pseudo_count_matrix,
                 threshold_fg, original_fg, original_bg, minimum_length, log_factorial_map, medium_pep_length, ZSCORE_FLAG):
        # All rp_pairs identified in this motif so far, NT with positive position, CT with negative position
        self.motif_rp_pairList = [] 
        self.num_rp = 0 # check if the rp is the first one of the motif, which will indicate the Stop sign
        self.MAX_LENGTH = max_length
        self.MINIMUM_LENGTH = minimum_length
        self.ZSCORE_FLAG = ZSCORE_FLAG
        self.BINOMIAL_PROB_THRESHOLD = binomial_prob_threshold
        self.NT_pseudo_count_matrix = NT_pseudo_count_matrix
        self.CT_pseudo_count_matrix = CT_pseudo_count_matrix
        self.threshold_fg = threshold_fg
        self.original_fg = original_fg
        self.original_bg = original_bg
        self.log_factorial_map = log_factorial_map
        self.deleterious_pairs_record = set([]) #keep record for deleterious pairs
        self.medium_pep_length = medium_pep_length
        
    # == Function Summary
    #   Perform iterative motif finding algorithm on input dataset and output
    #   one motif derived
    #   Loop over the two dataset until no significant peptide/position pairs could be identified any more 
    #   Two updated peptide datasets and one updated binomial log_prob matrix will be built
    # -- parameters 
    #    n_rp - number of recursion
    #    l_fold  - enrichment fold of current round of motif
    #    num_rp  - number of rp has been identified for the current motif
    #    pass_pairs - the pairs has been identified for the current motif
    def produce_one_motif(self, fg_branch, foreground_peps_set, bg_branch, background_peps_set, 
                          Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set):
        
        Stop = 0 # initilize the stop sign for the whole search
        
        sigPair_start_time = time.time() # Check run time of the program
        
        sigPair = fsp.SigPair(self.MAX_LENGTH, self.BINOMIAL_PROB_THRESHOLD, self.NT_pseudo_count_matrix,
                            self.CT_pseudo_count_matrix, self.motif_rp_pairList, alias_fg_set, 
                            alias_bg_set, self.log_factorial_map, self.medium_pep_length, self.ZSCORE_FLAG)
        
        self.num_rp += 1 # record the round of residue/position pair search
        [Flag, motif_aa, motif_pos, logp_sig_pair] = sigPair.find_one_pair(self.deleterious_pairs_record)
        sigPair_end_time = time.time()
        sigPair_runtime = (sigPair_end_time - sigPair_start_time)/60
        #print "sigPair time = ", sigPair_runtime, "min"
        
        # The number of fg and bg compared with the branch
        fg_vs_branch = (len(foreground_peps_set), len(fg_branch))
        bg_vs_branch = (len(background_peps_set), len(bg_branch))
        
        # Check if all significant pair has been found
        # and check if background is too small to have statistical significance
        if Flag == 1:
            if self.num_rp == 0 or not self.motif_rp_pairList:
                print("\nAll significant motifs have been identified and the search is exhaustive, done!")
                Stop = 1
            else:
                print("\nAll significant residue and position pairs have been found for the current motif")
            return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, 
                    last_p_sig_pair, l_fold, score, Stop)
        elif Flag == 2:
            print("Warning: Background become statistically unsatisfied, search stopped!!!")
            return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch,
                     last_p_sig_pair, l_fold, score, Stop)
        elif logp_sig_pair == 0:
            print ("Enforced stop the run because 0 binomial probability found, check the program!")
            Stop = 1
            return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, 
                    fg_vs_branch, bg_vs_branch, last_p_sig_pair, l_fold, score, Stop)
           
        # Keep tracking of recursion round
        n_rp += 1
            
        # ----------------------------------Iterative search and rebuild dataset and background set-------------
        # Iteratively search foreground dataset and rebuild both fore/background datasets by taking the 
        # regular expression from both datasets
            
        count_fg_find = 0 # count the number of peptide sequences removed from foreground peptide dataset
        
        # Initialize two sets to store updated fg and bg     
        new_foreground_set = set([])
        new_background_set = set([])
        new_alias_fg_set = set([])
        new_alias_bg_set = set([])
        
        # Store previous fg_branch and bg_branch for trace back one level
        backup_fg_branch = fg_branch
        backup_bg_branch = bg_branch
        
        # alias branch
        alias_fg_branch = alias_fg_set
        alias_bg_branch = alias_bg_set
        
        # Store previous bg and fg in datasets for fold enrichment
        fg_branch = foreground_peps_set
        bg_branch = background_peps_set
        old_fg_size = len(foreground_peps_set)
        old_bg_size = len(background_peps_set)
        
        for peptide in foreground_peps_set:
            pep_len = len(peptide)
            if (pep_len>motif_pos and motif_pos>=0):
                if peptide[motif_pos] == motif_aa:
                    new_foreground_set.add(peptide)
                    count_fg_find += 1
            elif motif_pos<0 and len(peptide)>=-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    new_foreground_set.add(peptide)
                    count_fg_find += 1
            else:
                continue  
            
        for peptide in alias_fg_set:
            pep_len = len(peptide)
            if (pep_len>motif_pos and motif_pos>=0):
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X' + peptide[motif_pos+1:]
                    new_alias_fg_set.add(alias_peptide)
            elif motif_pos == -1 and len(peptide)>-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X'
                    new_alias_fg_set.add(alias_peptide)
            elif motif_pos<-1 and len(peptide)>=-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X' + peptide[motif_pos+1:]
                    new_alias_fg_set.add(alias_peptide)
            else:
                continue       
        #update the foreground set by deducting the found motif patterns           
        foreground_peps_set = new_foreground_set
        # Number of sequences in the newest foreground peptide dataset
        fore_length = len(foreground_peps_set) 
        print "fg_set has ", fore_length, "peptides"
        alias_fg_set = new_alias_fg_set
        #print "alias_fg_set has ", len(alias_fg_set), "peptides"
        
        
                        
        for peptide in background_peps_set:
            pep_len = len(peptide)
            if (pep_len>motif_pos and motif_pos>=0):
                if peptide[motif_pos] == motif_aa:
                    new_background_set.add(peptide)
            elif motif_pos<0 and pep_len>=-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    new_background_set.add(peptide)
            else:
                continue
            
        for peptide in alias_bg_set:
            pep_len = len(peptide)
            if (pep_len>motif_pos and motif_pos>=0):
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X' + peptide[motif_pos+1:]
                    new_alias_bg_set.add(alias_peptide)
            elif motif_pos == -1 and len(peptide)>-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X'
                    new_alias_bg_set.add(alias_peptide)
            elif motif_pos<-1 and len(peptide)>=-motif_pos:
                if peptide[motif_pos] == motif_aa:
                    alias_peptide = peptide[:motif_pos] + 'X' + peptide[motif_pos+1:]
                    new_alias_bg_set.add(alias_peptide)
            else:
                continue         
            
        #update the background set by deducting the found motif patterns     
        background_peps_set = new_background_set
        # Number of sequences in the background peptide dataset
        back_length = len(background_peps_set) 
        print "bg_set has ", back_length, "peptides"
        alias_bg_set = new_alias_bg_set
        #print "alias_bg_set has ", len(alias_bg_set), "peptides"
		
        # keep record of fg and bg compared with the branch and use this when backtrack one-step of recursion
        kr_fg_vs_branch = fg_vs_branch
        kr_bg_vs_branch = bg_vs_branch
        
        # Update the number of fg and bg compared with the branch
        fg_vs_branch = (fore_length, len(fg_branch))
        bg_vs_branch = (back_length, len(bg_branch))
        
        # Calculate the Motif-X score, which equal to the sum of -log(binomial_p)
        score += -(logp_sig_pair) 

	#Check if identified rp pair is fault one, by adding the new position to existing rp position
        exception = 'False' # initialize var exception to indicate fault rp
        if self.motif_rp_pairList:
			if motif_pos >=0: # only check rp with opposite sign for conflict
				for rp in self.motif_rp_pairList:
                	#print "addition of two postions =", abs(motif_pos) + abs(rp[1])
					if rp[1] < 0 and abs(motif_pos) + abs(rp[1]) >= self.MAX_LENGTH:
						exception = 'True'
						break
			else:
				for rp in self.motif_rp_pairList:
                #print "addition of two postions =", abs(motif_pos) + abs(rp[1])
					if rp[1] >= 0 and abs(motif_pos) + abs(rp[1]) >= self.MAX_LENGTH:
						exception = 'True'
						break
        
        # each time the residue/position pair has been identified, mark it so won't count it in next round
        sig_pair = (motif_aa, motif_pos)
        self.motif_rp_pairList.append(sig_pair)
        # Check l_fold, i.e. the enrichment fold in the current cycle, to determine how to proceed
        
        # print out each recursion results for debugging
        print ("\n" +"top RP = "), motif_aa + ' / ' + str(motif_pos),
        print("(FG" + str(n_rp) + " = "), fore_length, 
        print(" , "),
        print("BG" + str(n_rp) + " = "), back_length,
        print(" ) ")
        print("Motif = "), self.motif_rp_pairList
        print (", binomial probability = "), exp(logp_sig_pair) # it will be zero if bp is too small
        last_p_sig_pair = logp_sig_pair
        #print('Z-Score = '), z_score
        #print(foreground_peps_set)
        
        l_fold = Motif.simple_fold_cal(self, fore_length, back_length, old_fg_size, old_bg_size) 
        print("Fold enriched current round = "), l_fold 
        
        # Check t_fold, i.e. the enrichment fold whole datawise
        # num_in_fg - total match of original peptides for the newest motif
        t_fold, num_in_fg, num_in_bg = utility.fold_enriched(self.original_fg, self.original_bg, self.motif_rp_pairList, "pairs")
        print("Fold enriched whole dataset = "), t_fold
        
        # If l_fold < 1 or 'NA', means deleterious fold or 0 background, 
        # if enrichment at current level less than 1, disregard the rp and back up one step            
        if l_fold <= 1 or t_fold <= 1 or exception=='True': 
            print("Warning: deleterious pair found, back up one step!")
            self.deleterious_pairs_record.add(sig_pair)
            #back up one step   
            (fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, 
                     fg_branch, bg_branch, score, n_rp, l_fold, alias_fg_set, alias_bg_set) = Motif.back_up_one_step(self, 
                                                                                        kr_fg_vs_branch,kr_bg_vs_branch, 
                        fg_branch, bg_branch, backup_fg_branch, backup_bg_branch, score, logp_sig_pair, n_rp, 
                        alias_fg_branch, alias_bg_branch)
            return Motif.produce_one_motif(self, fg_branch, foreground_peps_set, bg_branch, background_peps_set,
                                Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set)    
        
        elif fore_length < self.threshold_fg:
            if num_in_fg > self.threshold_fg:
                return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch,
                        logp_sig_pair, l_fold, score, Stop)
            else:
                print("Warining: foreground became too small to continue search")
                self.deleterious_pairs_record.add(sig_pair)
                # if foreground size become less than threshold, disregard the rp and back up one step
                (fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, 
                 fg_branch, bg_branch, score, n_rp, l_fold, alias_fg_set, alias_bg_set) = Motif.back_up_one_step(self, 
                        kr_fg_vs_branch,kr_bg_vs_branch, fg_branch, bg_branch, backup_fg_branch, backup_bg_branch, 
                        score, logp_sig_pair, n_rp, alias_fg_branch, alias_bg_branch) 
                
                # if the significant rp pair is too small, stop and return the current motif 
                return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch,
                        logp_sig_pair, l_fold, score, Stop)   
#                return Motif.produce_one_motif(self, fg_branch, foreground_peps_set, bg_branch, background_peps_set,
#                                Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set)   
        # motif won't take this rp, recursive call back to the next most significant rp 
        elif l_fold == 'NA':
            print('Background = 0, cannot calculate enrichment')
            if fore_length >= self.threshold_fg:
                return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch,
                        logp_sig_pair, l_fold, score, Stop)
                
            elif fore_length == 0:
                print "Motif positions collide, stop to debug!"
                Stop = 1
                return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch,
                        logp_sig_pair, l_fold, score, Stop)
            else:
                # if less than threshold or enrichment less than 2, disregard the rp and back up one step
#                if not Motif.significant_in_whole_set(self, t_fold, fore_length):
#                    print "Backup-one-step: motif not significant in whole dataset"
                self.deleterious_pairs_record.add(sig_pair)
                (fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, 
                    bg_branch, score, n_rp, l_fold, alias_fg_set, alias_bg_set) = Motif.back_up_one_step(self, 
                                                                                kr_fg_vs_branch,kr_bg_vs_branch,
                    fg_branch, bg_branch, backup_fg_branch, backup_bg_branch, score, logp_sig_pair, n_rp, 
                    alias_fg_branch, alias_bg_branch)
                
                return Motif.produce_one_motif(self, fg_branch, foreground_peps_set, bg_branch, background_peps_set,
                                Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set)
                         
        
        # If foreground_peps_set < set value of pseudo counts or background_peps_set < set value,
        # add pseudo counts to both matrix for calculation purpose
                 
#        elif (back_length < PSEUDO_N_SIZE / 2): 
#            print("Warning: Background less than half of " + str(PSEUDO_N_SIZE) + ", stop!")
#            return (self.motif_rp_pairList, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, 
#                    logp_sig_pair, l_fold, score, Stop)    
        
                    
                     
                 
        return Motif.produce_one_motif(self, fg_branch, foreground_peps_set, bg_branch, background_peps_set,
                                Flag, n_rp, last_p_sig_pair, l_fold, score, alias_fg_set, alias_bg_set)
        

    # == Function Summary
    #   Calculate the fold change for current motif by dividing the sequences matched a cetain regex motif
    #   in foreground set to the matched motif in background set
    #   Parameters:
    #   -- prompt, a string to give the correct file path
    #   Return:
    #   --fold enriched, if it is less than 1, indicating decreased fold
    def simple_fold_cal(self, new_fg_size, new_bg_size, old_fg_size, old_bg_size):
        fraction_fg = new_fg_size/float(old_fg_size)
        fraction_bg = new_bg_size/float(old_bg_size)
        if fraction_bg == 0:
            print("Warining: background is 0, could not calculate fold")
            fold = 'NA'
        else:
            fold = fraction_fg/fraction_bg
        return fold
          
    # == Function Summary
    #   Backtrack one step for current motif, the steps includes:
    #   1. Remove the significant residue/position from the current motif
    #   2. Backtrack one step to last fg_vs_branch, bg_vs_branch, fg, bg, fg_branch, bg_branch
    #   3. Backtrack to the last score by subtracting the current binomial possibility log value
    #   Parameters:
    #   -- kr_fg_vs_branch, kr_bg_vs_branch: record of last round of values
    #   -- fg_branch, bg_branch
    #   -- backup_fg_branch, backup_bg_branch
    #   -- score: score of motif
    #   -- logp_sig_pair: 
    #   -- n_rp: number of recursions
    #   Return:
    #   -- the group of values of last round 
    #      (pep_motif, fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, bg_branch, score, n_rp)
    #    
    def back_up_one_step(self, kr_fg_vs_branch, kr_bg_vs_branch, fg_branch, bg_branch, backup_fg_branch, backup_bg_branch,
                         score, logp_sig_pair, n_rp, alias_fg_branch, alias_bg_branch):
        
        self.motif_rp_pairList.pop()
        # backtrack one-step of recursion
        fg_vs_branch = kr_fg_vs_branch
        bg_vs_branch = kr_bg_vs_branch
        foreground_peps_set = fg_branch
        background_peps_set = bg_branch
        alias_fg_set = alias_fg_branch
        alias_bg_set = alias_bg_branch
        fg_branch = backup_fg_branch
        bg_branch = backup_bg_branch
        # Retro-calculate the Motif-X score, back up one step
        score -= -(logp_sig_pair)
        # The number of recursion round decrease by 1
        n_rp -= 1
        # Backtrack the enrichment fold
        l_fold = Motif.simple_fold_cal(self, len(foreground_peps_set), len(background_peps_set), 
                                       len(fg_branch), len(bg_branch)) 
        
        return (fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, 
                fg_branch, bg_branch, score, n_rp, l_fold, alias_fg_set, alias_bg_set)
    
    # == Function Summary
    #   Check if the current motif is significant in the whole dataset, 
    #   1. Match motif to the whole fg and bg dataset
    #   Parameters:
    #   -- pep_motif: motif with current rp
    #   -- original_fg, original_bg
    #   -- threshold_fg
    #   Return:
    #   -- True if it is significant else return false
    #    
    def significant_in_whole_set(self, t_fold, fore_length):
        if fore_length < self.threshold_fg: 
            print("Warning: Foreground with this motif has less than " + str(self.threshold_fg) + "in this round. ")
            print('Total peptides in foreground with this motif: '), 
            return False
        elif t_fold <= 1:
            print('Warning: Total enrichment = '), t_fold
            print("Iteration backtrack one level up and current residue/position pair disregarded due to fold enrichment")
            return False
        else:
            return True
            
        
