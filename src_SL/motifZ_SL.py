#!usr/bin/env python

# == Summary
#   Name: Motif-Z on same length peptides
#   New from version 1: add z-score to filter the significant residue-position pair

#   Date: August 20, 2013
#
# == Examples
#   This command generate peptide motifs with 9-mers:
#     python2.7 motif-Z-SL.py --inputf <filename> --ft [csv/csvGFY/txt] --outputf <filename> --length 9 
#
#   Like the first example, but with verbose/debug output:
#     python2.7 motif-Z-SL.py --inputf <filename> --ft [csv/csvGFY/txt] -V
#
#   This command prints out help info:
#     python2.7 motif-Z-SL.py -h
#
# == Usage
#   Usage: python2.7 motif-Z-SL.py --inputf <filename> --ft [csv/csvGFY/txt] --outputf <filename> --length 9 --bp 1e-7 
#
#   For help use: python2.7 motif-Z-SL.py -h
#
# == Options
#   -h, --help          Displays help message
#   -V, --verbose       Verbose output
#
# == Author
#   Xueheng Zhao


import re
import csv
import os
import sys
import copy
import time
from types import *
import timeit
from math import log
import utility
from decimal import *
from operator import itemgetter

from readFile_SL import *
import buildMatrix_SL as bmx
from utility import *
import buildMotif_SL as bm
import logging


# Global variables
global debug
debug = False



# == Function Summary
#   Main function to read in configuration file parameters, input datasets and output
#   motif results derived   
    
def runScript(file_type, peps, background_peps, result_file_path, debug_file_path, pep_length, zscore_flag,
              binomial_prob_threshold, threshold_fg):
    
    logging.basicConfig(filename='/var/www/motif-Z/log/motif_Z_SL.log', filemode='w', level=logging.DEBUG)
    logging.error('Started')
    start_time = time.time() # Check run time of the program
               
    """Background file are saved in .txt file with certain length (i.e. 7, 8, 9, 10, 11, 12, 20, 23) peptide sequences
    in each line by psedo-generated from peptide database"""
    
    #options = {7: '/var/www/motif-Z/bg_files/pep_bg_flank_7mers_100K.txt',
    #           8: '/var/www/motif-Z/bg_files/pep_bg_8mers_100K.txt',
    #           9: '/var/www/motif-Z/bg_files/pep_bg_9mers_1Million.txt',
    #           10: '/var/www/motif-Z/bg_files/pep_bg_10mers_100K.txt',
    #           11: '/var/www/motif-Z/bg_files/pep_bg_11mers_100K.txt',
    #           12: '/var/www/motif-Z/bg_files/pep_bg_12mers_100K.txt',
    #           20: '/var/www/motif-Z/bg_files/pep_bg_20mers_100K.txt',
    #           23: '/var/www/motif-Z/bg_files/pep_bg_23mer_100K_wFlank.txt'
    #}
    
    #background_file = options[int(pep_length)]

    # Constants
    MOTIF_LENGTH = int(pep_length)
    AA_LIST = 'ARNDCEQGHILKMFPSTWYV'
    Z_LIST = ['1', '2', '3', '4', '5', '6']
    
    if zscore_flag in Z_LIST:
        zscore = int(zscore_flag)
        ZSCORE_FLAG = zscore
        BINOMIAL_PROB_THRESHOLD = log(float(1000))
    else:
        ZSCORE_FLAG = 0
        BINOMIAL_PROB_THRESHOLD = log(float(binomial_prob_threshold))
    #BINOMIAL_PROB_THRESHOLD = float(binomial_prob_threshold)
    #print('BINOMIAL_PROB_THRESHOLD = '), BINOMIAL_PROB_THRESHOLD
    
    
    Flag = 0
    PSEUDO_N_SIZE = 50      
        
# == Data structure used to store FG and BG peptide sequences: Set
    
    # Foreground peptide array save in the data structure set
    foreground_peps_set = set([])
    num_fg_pep = 0 # number of peptide seqs in the foreground
    for pep in peps:
        if not '.' in pep:
            peptide = re.sub("[^A-Za-z]", "", pep)
            if len(peptide)== MOTIF_LENGTH:
                    num_fg_pep += 1
                    foreground_peps_set.add(peptide)
    peps_size = len(foreground_peps_set)
    print 'Input files contains peptides with length ', MOTIF_LENGTH, ' = ', peps_size
    
        
    # Build background dataset from input peptide file
    background_peps_set = set([])  # Background tpeptide array (control)
    #num_bg_pep = 0 # number of motifs in the background
    for b_pep in background_peps:
        b_pep = b_pep.upper()
        b_peptide = re.sub("[^ARNDCEQGHILKMFPSTWYV]", "U", b_pep)
        if len(b_peptide) == MOTIF_LENGTH:
     #       num_bg_pep += 1
            background_peps_set.add(b_peptide)


    # Check the runtime of data readin
    #print('runtime of reading data = '), (time.time() - start_time)/60.0, "min"
    start_time = time.time() 
    
    # Save a copy of both fore/background datasets since the processing carried out in this step modify the
    # original dataset lists   
        
    original_fg = copy.deepcopy(foreground_peps_set)
    original_bg = copy.deepcopy(background_peps_set)
 
    
    # Size of foreground and background 
    fg_size = len(original_fg)
    bg_size = len(original_bg)
    #print('fg_size = '), fg_size
    #print('bg_size = '), bg_size
    print 'current run in progress....'
    
    # Calculate pseudo count matrix
    pc_matrix = bmx.Matrix(MOTIF_LENGTH, original_bg)
    pseudo_count_matrix = pc_matrix.aaFreqMatrix() 
    # Build the precomputated factorial map
    log_factorial_map = utility.log_factorial_map(fg_size)
    
       

    # Set up a data structure for table output about size in fg and bg in previous motif, and binomial probability of last 
    # residue/position pair
    motif_pre_ref = {} 
    r_count = 0 # Initialize the recursion counts, this number added by itself in the recursive call
    all_motif = set([]) # Initialize motif set to store the motifs identified from input data
    
    # Check the runtime of data readin
    #print('runtime of build one set of matrices = '), (time.time() - start_time)/60.0, "min"
    start_time = time.time()
    
    # Save stand output into the debug file with format for debugging
    saveout = sys.stdout
    fhout = open(debug_file_path, "w")      
    sys.stdout = fhout 
    
    
    print("---------------------------------------------------------------------------------")
    print("Start recursion")
    
    all_motifs_found, motif_pre_ref = bm.produce_all_motif(MOTIF_LENGTH, foreground_peps_set, background_peps_set,
                BINOMIAL_PROB_THRESHOLD, Flag, all_motif, r_count, pseudo_count_matrix, motif_pre_ref, threshold_fg, 
                original_fg, original_bg, log_factorial_map, ZSCORE_FLAG)
    
    # Save stand output to debug file and close the file handle
    sys.stdout = saveout                                     
    fhout.close()
    
    # Check the runtime of data readin
    #print('runtime of recursive search = '), (time.time() - start_time)/60.0, "min"
    start_time = time.time()
    
    # Write results in the result file
    f_result = open(result_file_path, 'wb')  
    f_result.write("Foreground has total peptides = " + str(fg_size) + "\n") 
    f_result.write("Backgorund has total peptides = " + str(bg_size) + "\n") 
     
    if os.path.exists(result_file_path):
        result_data = {} # The motifs identified and the fold change results are stored in a dictionary structure
        for motif in all_motifs_found:
            fold_enrichment, num_in_fg, num_in_bg = bm.fold_enriched(original_fg, original_bg, motif)
            if fold_enrichment > 1:
                result_data[motif] = (fold_enrichment, num_in_fg, num_in_bg)

        
        # write the result in file with str.format() 
        f_result.write('-----------------------------------------' + '---------------------------------------' + 
                       '-----------------------------------------' + '--------------------------------' + '\n') 

        f_result.write('{0:25}\t {1:7}\t {2:7}\t {3:7}\t {4:7}\t {5:7}\t {6:8}\t {7:8}\t {8:9}\n'.format("motif", "t_enrichment", "in_t_fg", 
                        "fg_vs_branch", "in_t_bg", "bg_vs_branch", "log_p", "l_enrichment", "score"))
        f_result.write('-----------------------------------------' + '---------------------------------------' + 
                       '-----------------------------------------' + '--------------------------------' + '\n')   
        for motif in result_data:
            # Format the bg_vs_branch and fg_vs_branch output
            motif_pre_ref[motif][0] = str(motif_pre_ref[motif][0][0]) + '/' + str(motif_pre_ref[motif][0][1])
            motif_pre_ref[motif][1] = str(motif_pre_ref[motif][1][0]) + '/' + str(motif_pre_ref[motif][1][1])
            # Wirte output data into the format described below
            f_result.write('{0:25}\t {1:7}\t {2:7}\t {3:7}\t {4:7}\t {5:7}\t {6:8.2e}\t {7:7.1f}\t {8:8.1f}\n'.format(motif, 
                            result_data[motif][0], result_data[motif][1], motif_pre_ref[motif][0], result_data[motif][2], 
                            motif_pre_ref[motif][1], motif_pre_ref[motif][2], motif_pre_ref[motif][3], motif_pre_ref[motif][4]))
        f_result.close()
    else:
        logging.error("The file " + result_file_path + " doesn't exist")
        
    # Build a dictionary to store motif regular expression results and use this dictionary in the webpage output
    result_table = []
    for rexp_motif in result_data:
        total_enrich = result_data[rexp_motif][0]
        match_fg = result_data[rexp_motif][1]
        p_value = '%.2E' % exp(motif_pre_ref[rexp_motif][5])
        score = motif_pre_ref[rexp_motif][4]
	rpSeqList = rp_iteration_id(motif_pre_ref[rexp_motif][6])
        result_table.append([rexp_motif, rpSeqList, total_enrich, match_fg, p_value, score])
    # Sort the table based on score
    sorted_result_table = sorted(result_table, key=itemgetter(5), reverse=True)
    # Format score with 2 dicemal place
    for item in sorted_result_table:
	item[5] = '%.2f' % item[5]

    logging.error('Finished') 
    return sorted_result_table
