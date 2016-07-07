    # == Summary
#   This module contains the major iteration and recursion procedure for Motif-Z-SL algorithm
#   Perform iterative motif finding algorithm on input file and writes output
#   
#
# == Functions
#  


import re
import buildMatrix_SL as bm
import utility
from math import exp
import checkZScore_SL as cz

# Global variables
global debug
debug = False

# Constant and default values
PSEUDO_N_SIZE = 50
AA_LIST = 'ARNDCEQGHILKMFPSTWYV'

# == Function Summary
#   One-iteration algorithm to find one significant residue/position pair
#   build the peptide motif from the iterative search
#   Parameters: 
#   -- foreground_peps_set, peptide dataset as foreground 
#   -- background_peps_set, peptide dataset as background
#   -- MOTIF_LENGTH, the length of peptide
#   -- BINOMIAL_PROB_THRESHOLD, value used as threshold to identify significant pair
#   -- pseudo_count_matrix, the probability matrix generated from the beginning background set
#   -- num_rp  - number of rp has been identified for the current motif
#   -- pass_pairs - the pairs has been identified for the current motif
#   Return:
#   function return a tuple of values including Flag, motif_aa, motif_pos, logp_sig_pair, fg_vs_branch, bg_vs_branch
#   example: (1, Y, 8, 1e-10, (25,100), (50,200))
#   
def find_one_pair(foreground_peps_set, background_peps_set, BINOMIAL_PROB_THRESHOLD, MOTIF_LENGTH, 
                  pseudo_count_matrix, threshold_fg, num_rp, fixed_pairs, 
                  log_factorial_map, deleterious_pairs_record, ZSCORE_FLAG):    
    
    Flag =0 # initialize flag to indicate the status of algorithm
    
    
    #----------------------------------Step 1, build matrices --------------------------------------
    # Build the dataset_frequency_matrix as a database using dictionary type
    fg_matrix = bm.Matrix(MOTIF_LENGTH, foreground_peps_set) 
    fg_count_matrix = fg_matrix.aaCountMatrix()
    # Build the background_probablity_matrix as a database using dictionary type
    bg_freq_matrix = (bm.bgMatrixWithPseudoCount(MOTIF_LENGTH, background_peps_set, 
                                                             pseudo_count_matrix).getMatrixWithPC())
        
    if len(background_peps_set) == 0:
        Flag = 2
        print("Background does not contain this motif, exit from search")
        return (Flag, '137', 137, 137, num_rp)
            
    # Calculate binomial matrix                                                      
    bpl_matrix = fg_matrix.logbpMatrix(fg_count_matrix, bg_freq_matrix, log_factorial_map)
    
    # Generate the significant rp pair list based on logbp values
    sig_pair_list = sigPairList(bpl_matrix, MOTIF_LENGTH, BINOMIAL_PROB_THRESHOLD, fixed_pairs)
    if len(sig_pair_list) == 0:
        # No more significant residue/position pairs are discovered
        Flag = 1
        return (Flag, '9', 137, 0, num_rp)
    
    test_aa, test_pos, test_logbp = findSigPair(sig_pair_list, fg_count_matrix, 
                                                deleterious_pairs_record, MOTIF_LENGTH)
    if test_aa == 'empty':
        Flag = 1
        return (Flag, '9', 137, 0, num_rp)        
    
    print "sig_pair = ", test_aa, test_pos, test_logbp
    
    if ZSCORE_FLAG == 0:
        return (Flag, test_aa, test_pos, test_logbp, num_rp)
    else:
        checkZ = cz.CheckZScore(MOTIF_LENGTH, foreground_peps_set, fixed_pairs, bg_freq_matrix,
                                     test_logbp, log_factorial_map) 
        z_score = checkZ.calculate_ZScore()
    
        print("z-score = "), z_score
    
        if z_score > ZSCORE_FLAG:
            # print out the significant bp value for debug file  
            rp_logbp = test_logbp 
            rp_aa = test_aa
            rp_pos = test_pos
            print("log form of binomial probability of significant rp pair = "), rp_logbp
        else:
            print("This pair doesn't pass z_score test, return this motif")
            Flag = 1
            return (Flag, '9', 137, 0, num_rp)
        
    num_rp += 1 # record the round of residue/position pair search
                
    return (Flag, rp_aa, rp_pos, rp_logbp, num_rp)
    
    
# Traverse the binomial_log_matrix to determine the best hits based on the probability less 
# than the threshold. Data structure used is a dictionary with a tuple (aa, pos) as key and 
# the binomial probability as value
def sigPairList(bplog_matrix, MOTIF_LENGTH, BINOMIAL_PROB_THRESHOLD, fixed_pairs):   
    sig_pair_list = []
    for aa in AA_LIST:
        for pos in range(MOTIF_LENGTH):
            sig_pair_dict = {}
            # Check if the binomial probability == 0, this caused by either p == 0, which is background have 0 r/p pair
            # or p == 1, which is the case after pruning
            if (bplog_matrix[aa][pos] < BINOMIAL_PROB_THRESHOLD and bplog_matrix[aa][pos] != 0):
                sig_pair = (aa, pos)
                if sig_pair not in fixed_pairs:
                    sig_pair_dict['pair'] = sig_pair
                    sig_pair_dict['log_prob'] = bplog_matrix[aa][pos]
                    sig_pair_list.append(sig_pair_dict)
                        
    return sig_pair_list
    
# Find significant pair in one list if another is empty
def findSigPair(sig_pair_list, fg_count_matrix, deleterious_pairs_record, MOTIF_LENGTH):
    
    #Check the pairs have been tried
    pairs_to_remove = list()
        
    rp_aa, rp_pos, rp_logbp = ('debug', 'debug', 'debug')
        
    for rp_dict in sig_pair_list:
        if rp_dict['pair'] in deleterious_pairs_record:
            print "rp_dict['pair'] removed= ", rp_dict['pair']
            pairs_to_remove.append(rp_dict)
    if len(pairs_to_remove) > 0:
        for rp_dict in pairs_to_remove:
            sig_pair_list.remove(rp_dict)
            
    print "deleterious_pairs_record=", deleterious_pairs_record
                    
    # check if either list is empty
    if len(sig_pair_list) == 0:
        return 'empty', 'empty', 'empty'
    elif len(sig_pair_list) == 1:
    # The only one residue and position pair having significant log_prob
        rp_aa = sig_pair_list[0]['pair'][0]
        rp_pos = sig_pair_list[0]['pair'][1]
        rp_logbp = sig_pair_list[0]['log_prob']
    else:         
        # Find residue aa and position pair(s) having the lowest binomial log_prob            
        sig_pair_list.sort(lambda x, y : cmp(x['log_prob'], y['log_prob']))
        sorted_sig_pair_list = sig_pair_list
        #print "sorted_sig_pair_list=", sorted_sig_pair_list
              
        if sorted_sig_pair_list[0]['log_prob'] < sorted_sig_pair_list[1]['log_prob']:
            rp_aa = sorted_sig_pair_list[0]['pair'][0]
            rp_pos = sorted_sig_pair_list[0]['pair'][1]       
            rp_logbp = sorted_sig_pair_list[0]['log_prob']
        else:
            # list to save the pairs with equal binomial log_prob
            motif_pairs = [sorted_sig_pair_list[0]['pair']]  
            for i in range(1, len(sorted_sig_pair_list)):
                if sorted_sig_pair_list[i]['log_prob'] == sorted_sig_pair_list[0]['log_prob']:
                    motif_pairs.append(sorted_sig_pair_list[i]['pair'])
                        
            # chose the aa with the lowest log_prob and highest count in the dataset frequency matrix 
            mr = 0 
                
            for j in range(len(motif_pairs)):
                if (fg_count_matrix[motif_pairs[j][0]][motif_pairs[j][1]] > 
                fg_count_matrix[motif_pairs[mr][0]][motif_pairs[mr][1]]):
                    mr = j
            rp_aa = motif_pairs[mr][0]
            rp_pos = motif_pairs[mr][1]        
            rp_logbp = sorted_sig_pair_list[0]['log_prob']
            
    return rp_aa, rp_pos, rp_logbp

        
# == Function Summary
#   Perform iterative motif finding algorithm on input dataset and output
#   one motif derived
# -- parameters 
#    n_rp - number of recursion
#    l_fold  - enrichment fold of current round of motif
#    num_rp  - number of rp has been identified for the current motif
#    pass_pairs - the pairs has been identified for the current motif
def produce_one_motif(MOTIF_LENGTH, fg_branch, foreground_peps_set, bg_branch, background_peps_set, BINOMIAL_PROB_THRESHOLD, 
                      Flag, pep_motif, n_rp, pseudo_count_matrix, last_p_sig_pair, l_fold, score, threshold_fg, 
                      original_fg, original_bg, num_rp, log_factorial_map, deleterious_pairs_record, fixed_pairs, ZSCORE_FLAG):
    
    # Loop over the two dataset until no significant peptide/position pairs could be identified any more 
    # Two updated peptide datasets and one updated binomial log_prob matrix will be built
    
    Stop = 0 # initilize the stop sign for the whole search
    [Flag, motif_aa, motif_pos, logp_sig_pair, num_rp] = find_one_pair(foreground_peps_set, background_peps_set, 
                                                                                BINOMIAL_PROB_THRESHOLD, MOTIF_LENGTH, 
                                                                                pseudo_count_matrix, threshold_fg, num_rp, 
                                                                                fixed_pairs, 
                                                                                log_factorial_map, deleterious_pairs_record, 
                                                                                ZSCORE_FLAG)
    
    # The number of fg and bg compared with the branch
    fg_vs_branch = (len(foreground_peps_set), len(fg_branch))
    bg_vs_branch = (len(background_peps_set), len(bg_branch))
    
    # Check if all significant pair has been found
    # and check if background is too small to have statistical significance
    if Flag == 1:
        if num_rp == 0:
            print("\nAll significant motifs have been identified and the search is exhaustive, done!")
            Stop = 1
        else:
            print("\nAll significant residue and position pairs have been found for the current motif")
        return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, last_p_sig_pair, l_fold, score, Stop, fixed_pairs)
    elif Flag == 2:
        print("Warning: Background become statistically unsatisfied, search stopped!!!")
        return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, last_p_sig_pair, l_fold, score, Stop, fixed_pairs)
    elif logp_sig_pair == 0:
        print ("Enforced stop the run because 0 binomial probability found, check the program!")
        Stop = 1
        for i in range(MOTIF_LENGTH):
            pep_motif = pep_motif.replace('.', 'X')
        return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, last_p_sig_pair, l_fold, score, Stop, fixed_pairs)
    
    # Check z-Score before adding to the motif
    
    #update the motif
    pep_motif = pep_motif[:motif_pos] + motif_aa + pep_motif[motif_pos+1:]
        
    # Keep tracking of recursion round
    n_rp += 1
      
    #------------------------------------Regular Expression------------------------------------------ 
    # Build the motif in regex format
        
    # Compile the motif into regex format
    regex_motif = re.compile(pep_motif)
        
    # ----------------------------------Iterative search and rebuild dataset and background set-------------
    # Iteratively search foreground dataset and rebuild both fore/background datasets by taking the 
    # regular expression from both datasets
        
    count_fg_find = 0 # count the number of peptide sequences removed from foreground peptide dataset
    
    # Initialize two sets to store updated fg and bg     
    new_foreground_set = set([])
    new_background_set = set([])
    
    # Store previous fg_branch and bg_branch for trace back one level
    backup_fg_branch = fg_branch
    backup_bg_branch = bg_branch
    
    # Store previous bg and fg in datasets for fold enrichment
    fg_branch = foreground_peps_set
    bg_branch = background_peps_set
    old_fg_size = len(foreground_peps_set)
    old_bg_size = len(background_peps_set)
    
    for peptide in foreground_peps_set:
        if regex_motif.match(peptide):
            new_foreground_set.add(peptide)
            count_fg_find += 1 
    #update the foreground set by deducting the found motif patterns           
    foreground_peps_set = new_foreground_set
    fore_length = len(foreground_peps_set) # Number of sequences in the foreground peptide dataset
    
                    
    for peptide in background_peps_set:
        if regex_motif.match(peptide):
            new_background_set.add(peptide)  
    #update the background set by deducting the found motif patterns     
    background_peps_set = new_background_set
    back_length = len(background_peps_set) # Number of sequences in the background peptide dataset
    
    # keep record of fg and bg compared with the branch and use this when backtrack one-step of recursion
    kr_fg_vs_branch = fg_vs_branch
    kr_bg_vs_branch = bg_vs_branch
    
    # Update the number of fg and bg compared with the branch
    fg_vs_branch = (len(foreground_peps_set), len(fg_branch))
    bg_vs_branch = (len(background_peps_set), len(bg_branch))
    
    # print out each recursion results for debugging
    print ("\n" +"top RP = "), motif_aa + ' / ' + str(motif_pos),
    print("(FG" + str(n_rp) + " = "), len(foreground_peps_set), 
    print(" , "),
    print("BG" + str(n_rp) + " = "), len(background_peps_set),
    print(" ) ")
    print("Motif = "), pep_motif,
    print (", binomial probability = "), exp(logp_sig_pair)
    last_p_sig_pair = logp_sig_pair
    #print('Z-Score = '), z_score
    #print(foreground_peps_set)
    
    # Calculate the Motif-X score, which equal to the sum of -log(binomial_p)
    score += -(logp_sig_pair) 
    
    # each time the residue/position pair has been identified, mark it so won't count it in next round
    pass_pair = (motif_aa, motif_pos)
    fixed_pairs.append(pass_pair)
    
    # Check l_fold, i.e. the enrichment fold in the current cycle, to determine how to proceed
    l_fold = simple_fold_cal(fore_length, back_length, old_fg_size, old_bg_size) 
    print("Fold enriched current round = "), l_fold 
    
    # Check t_fold, i.e. the enrichment fold whole datawise
    t_fold, num_in_fg, num_in_bg = utility.fold_enriched(original_fg, original_bg, regex_motif, "rexp")
    print("Fold enriched whole dataset = "), t_fold

    
    # If l_fold < 1 or 'NA', means deleterious fold or 0 background, 
    # motif won't take this rp, recursive call back to the next most significant rp 
    if l_fold == 'NA':
        print('Background = 0, cannot calculate enrichment')
        if fore_length >= threshold_fg or significant_in_whole_set(pep_motif, original_fg, original_bg, threshold_fg):
            return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, logp_sig_pair, l_fold, 
                    score, Stop, fixed_pairs)
        else:
            # if less than threshold or enrichment less than 2, disregard the rp and back up one step
            (pep_motif, fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, bg_branch, score, 
                 n_rp, l_fold) = back_up_one_step(pep_motif, motif_pos, kr_fg_vs_branch,kr_bg_vs_branch, fg_branch, bg_branch, backup_fg_branch, 
                                          backup_bg_branch, score, logp_sig_pair, n_rp)
            return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, logp_sig_pair, l_fold, 
                        score, Stop, fixed_pairs)
                     
            
    # if enrichment at current level less than 1, disregard the rp and back up one step            
    elif l_fold <= 1 or t_fold <= 1: 
        
        print("Warning: deleterious pair found, fold enrichment = " + str("{0:.2}".format(l_fold)))
        # back up one step   
        (pep_motif, fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, bg_branch, score, 
                 n_rp, l_fold) = back_up_one_step(pep_motif, motif_pos, kr_fg_vs_branch,kr_bg_vs_branch, fg_branch, bg_branch, 
                                                  backup_fg_branch, backup_bg_branch, score, logp_sig_pair, n_rp)
        fixed_pairs.pop()
 
    # If foreground_peps_set < set value of pseudo counts or background_peps_set < set value,
    # add pseudo counts to both matrix for calculation purpose
        
    elif fore_length < threshold_fg:
        # if in this motif has more than threshold in fg, keep it
        if num_in_fg >= threshold_fg:
            return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, logp_sig_pair, l_fold, 
                score, Stop, fixed_pairs)
        else:
            # if matched peptides in original foreground less than threshold, disregard the rp and back up one step
            (pep_motif, fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, bg_branch, score, 
                 n_rp, l_fold) = back_up_one_step(pep_motif, motif_pos, kr_fg_vs_branch,kr_bg_vs_branch, fg_branch, bg_branch, backup_fg_branch, 
                                          backup_bg_branch, score, logp_sig_pair, n_rp)
	    fixed_pairs.pop()
            return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, logp_sig_pair, l_fold, 
                score, Stop, fixed_pairs)
    
    elif (len(background_peps_set) < PSEUDO_N_SIZE / 2): 
        print("Warning: Background less than half of " + str(PSEUDO_N_SIZE) + ", stop!")
        return (pep_motif, foreground_peps_set, background_peps_set, n_rp, fg_vs_branch, bg_vs_branch, logp_sig_pair, l_fold, 
                score, Stop, fixed_pairs)
            
    elif fore_length / back_length > 2:
        print("Warning: Background less than half of foreground!") 
                
    return produce_one_motif(MOTIF_LENGTH, fg_branch, foreground_peps_set, bg_branch, background_peps_set, BINOMIAL_PROB_THRESHOLD, 
                                         Flag, pep_motif, n_rp, pseudo_count_matrix, last_p_sig_pair, l_fold, score, threshold_fg,
                                         original_fg, original_bg, num_rp, log_factorial_map, deleterious_pairs_record, fixed_pairs, ZSCORE_FLAG)
    
# == Function Summary
#   Perform iterative motif finding algorithm on input dataset and output motif set derived
#   recursion will update both foreground and background based on each motif identified in previous search circle 

def produce_all_motif(MOTIF_LENGTH, foreground_peps_set, background_peps_set,BINOMIAL_PROB_THRESHOLD, 
                      Flag, all_motif, r_count, pseudo_count_matrix, motif_pre_ref, threshold_fg,
                       original_fg, original_bg, log_factorial_map, ZSCORE_FLAG):
    
    
    if r_count == 0:
        print("####################################################################################")
        print("Round " + str(r_count) + ": "),
        print("FG" + str(r_count) +" = "+ str(len(foreground_peps_set)) + "\t BG" + str(r_count) +
           " = " + str(len(background_peps_set)))
    
    # for debug purpose, stop at 10 recursion steps
    """if r_count == 10:
        print("for debugging, stop here")
        return all_motif"""
    
    
    # Loop over until no peptide motif can be identified  
    """if Flag == 1:
        print("Note: All significant residue and position pairs have been found for the input data file")
        return all_motif, motif_pre_ref"""
    if len(background_peps_set) == 0:
        print("Warning: Empty background, suspend search!")
        return all_motif, motif_pre_ref
    elif len(foreground_peps_set) == 0:
        print("Warning: Empty foreground, suspend search!")
        return all_motif, motif_pre_ref
    
    # Initialize the motif using regex like format
    pep_motif = '' 
    for i in range(MOTIF_LENGTH):
        pep_motif += '.'
        
    # Make the checkpoint for all the significant r/p has exhausted
    blank = pep_motif
    
    # This variables just for debugging
    #bg0 = len(background_peps_set)
    start_fg = foreground_peps_set
    start_bg = background_peps_set
    
    # initialize fg_branch and bg_branch
    fg_branch = start_fg
    bg_branch = start_bg
    
    
    # Produce one motif at a time
    # This used for debugging
    n_rp = 0 # counting the number of residue/position pair found
    num_rp = 0 # check if the rp is the first one of the motif, which will indicate the Stop sign
    # parameter for output of binomial probability for the last residue/position pair, initialize to 1
    last_p_sig_pair = 1
    # initialize the score for each motif
    score = 0
    l_fold = 0
    # initialize deleterious pairs for tracking
    deleterious_pairs_record = set([])
    fixed_pairs=[]
            
    [motif_identified, searched_fg_set, searched_bg_set, debug_n, fg_vs_branch, bg_vs_branch, binomial_p, 
     l_enrichment, score, Stop, rpSeqList] = produce_one_motif(MOTIF_LENGTH, fg_branch, foreground_peps_set, 
            bg_branch, background_peps_set, BINOMIAL_PROB_THRESHOLD, Flag, pep_motif, n_rp, pseudo_count_matrix, 
        last_p_sig_pair, l_fold, score, threshold_fg, original_fg, original_bg, num_rp, log_factorial_map, 
        deleterious_pairs_record, fixed_pairs, ZSCORE_FLAG)
     
    # Check if all motifs have been exhausted, if so, stop recursion    
    if blank == motif_identified: 
        print("Note: All significant pairs identified and search is exhaustive for the input data file")  
        return all_motif, motif_pre_ref 

    # print out information for the debug file to track the fg and bg size
    print("Final Motif = "),
    print(motif_identified),
    print("rp pairs= "),
    print(rpSeqList),
    print("(" + str(len(searched_fg_set)) + "/ " + str(len(searched_bg_set))) + ") \t",
    motif_enrichment, num_in_fg, num_in_bg = fold_enriched(start_fg, start_bg, motif_identified)
    print("Fold enrichment = "),
    #Calculate p value for the motif, using binomial probability for the whole motif
    n = len(start_fg)
    k = num_in_fg
    p = float(num_in_bg)/len(start_bg)
    log_bp_motif = utility.binomial_log_prob(n, k, p, log_factorial_map)
    
    if motif_enrichment == 'NA':
        print('NA')
    else:
        print("{0:.2f}".format(motif_enrichment))
    
    #update the motif
    if motif_identified not in all_motif and motif_identified != blank:
        all_motif.add(motif_identified)
        
    # Save the reference value for current motif
    motif_pre_ref[motif_identified] = [fg_vs_branch, bg_vs_branch, binomial_p, l_enrichment, score, log_bp_motif, rpSeqList]
        
    
    print("Total Motif = "),
    print(all_motif),
    
    print("\t Total Motif size = "),
    print(len(all_motif))
    
    #Save the current foreground and background set as fg_branch and bg_branch
    fg_branch = foreground_peps_set
    bg_branch = background_peps_set
    
    #update the foreground and background sets by deducting the found motif patterns
    foreground_peps_set = foreground_peps_set - searched_fg_set
    background_peps_set = background_peps_set - searched_bg_set
    
    r_count += 1 # Recursion counting
    print("#################################################################################################")
    print("Round " + str(r_count) + ": "),
    print("FG" + str(r_count) + " = FG" + str(r_count - 1) + " - " + " FG" + str(debug_n) + " = " + str(len(foreground_peps_set)) + "\t BG" + str(r_count) +
           " = " + str(len(background_peps_set)))
    

    return produce_all_motif(MOTIF_LENGTH, foreground_peps_set, background_peps_set,
                BINOMIAL_PROB_THRESHOLD, Flag, all_motif, r_count, pseudo_count_matrix, motif_pre_ref, 
                threshold_fg, original_fg, original_bg, log_factorial_map, ZSCORE_FLAG)    
    
    
    
# == Function Summary
#   Count how many peptides in a dataset matched a certain regex motif
#   output the number of peptides in the dataset by the motif pattern match

def number_matched_motif(this_set, motif):
    
    count = 0 # initialize the number of motifs matched in the dataset
    regex_motif = re.compile(motif) # Compile the motif into regex format
    for item in this_set:
        if regex_motif.match(item):
            count += 1
    return count

# == Function Summary
#   Calculate the fold change for current motif by dividing the sequences matched a cetain regex motif
#   in foreground set to the matched motif in background set
#   Parameters:
#   -- prompt, a string to give the correct file path
#   Return:
#   --the value of the fold of enrichment

def fold_enriched(set1, set2, motif): 
    
    num_in_fg = number_matched_motif(set1, motif)
    num_in_bg = number_matched_motif(set2, motif)
    
    fg_size = len(set1)
    bg_size = len(set2)
    
    fraction_motif_set1 = num_in_fg / float(fg_size)
    fraction_motif_set2 = num_in_bg / float(bg_size)
    
    if fraction_motif_set2 != 0:
        fold = fraction_motif_set1 / fraction_motif_set2
        if fold > 2:
            fold = int(fold) 
        else:
            fold = round(fold, 2)  
    else:
        print("Warning: zero background in fold enrichment calculation")  
        fold = 'NA'
    return fold, num_in_fg, num_in_bg


# == Function Summary
#   Calculate the fold change for current motif by dividing the sequences matched a cetain regex motif
#   in foreground set to the matched motif in background set
#   Parameters:
#   -- prompt, a string to give the correct file path
#   Return:
#   --fold enriched, if it is less than 1, indicating decreased fold
def simple_fold_cal(new_fg_size, new_bg_size, old_fg_size, old_bg_size):
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
#   -- pep_motif: motif with current rp
#   -- motif_pos: position of rp
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
def back_up_one_step(pep_motif, motif_pos, kr_fg_vs_branch, kr_bg_vs_branch, fg_branch, bg_branch, backup_fg_branch, backup_bg_branch,
                     score, logp_sig_pair, n_rp):
    
    pep_motif = pep_motif[:motif_pos] + '.' + pep_motif[motif_pos+1:]
    # backtrack one-step of recursion
    fg_vs_branch = kr_fg_vs_branch
    bg_vs_branch = kr_bg_vs_branch
    foreground_peps_set = fg_branch
    background_peps_set = bg_branch
    fg_branch = backup_fg_branch
    bg_branch = backup_bg_branch
    # Retro-calculate the Motif-X score, back up one step
    score -= -(logp_sig_pair)
    # The number of recursion round decrease by 1
    n_rp -= 1
    # Backtrack the enrichment fold
    l_fold = simple_fold_cal(len(foreground_peps_set), len(background_peps_set), len(fg_branch), len(bg_branch)) 
    
    return (pep_motif, fg_vs_branch, bg_vs_branch, foreground_peps_set, background_peps_set, fg_branch, bg_branch, score, n_rp, l_fold)

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
def significant_in_whole_set(pep_motif, original_fg, original_bg, threshold_fg):
    num_t_match = number_matched_motif(original_fg, pep_motif)
    num_t_enrich = fold_enriched(original_fg, original_bg, pep_motif)
    if num_t_match < threshold_fg or num_t_enrich[0] < 2: 
        print("Warning: Foreground with this motif has less than " + str(threshold_fg) + "in this round. ")
        print('Total peptides in foreground with this motif: '), num_t_match
        print('Total enrichment = '), num_t_enrich[0]
        print("Iteration backtrack one level up and current residue/position pair disregarded")
        return False
    else:
        return True
        
        
        
