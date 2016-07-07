'''
Modified on July 30, 2013

@author: zhaoxh
'''

import buildMatrix as bm  
import utility
from math import exp
import time
import checkZScore

# Global variables
global debug
debug = False

# Constant and default values
PSEUDO_N_SIZE = 50
AA_LIST = 'ARNDCEQGHILKMFPSTWYV'

class SigPair(object):
    '''
    One-iteration algorithm to find one significant residue/position pair
    '''
    def __init__(self, max_length, binomial_prob_threshold, NT_pseudo_count_matrix, CT_pseudo_count_matrix,
                 iPairSet, fg_data_set, bg_data_set, log_factorial_map, medium_pep_length, Z_flag):
        '''
        Constructor for instance parameters
        '''
        self.mypair = ()
        self.iPairSet = iPairSet
        self.MAX_LENGTH = max_length
        self.BINOMIAL_PROB_THRESHOLD = binomial_prob_threshold
        self.Z_flag = Z_flag
        self.NT_pseudo_count_matrix = NT_pseudo_count_matrix
        self.CT_pseudo_count_matrix = CT_pseudo_count_matrix
        self.fg_data_set = fg_data_set # alias_fg_set
        self.bg_data_set = bg_data_set # alias_bg_set
        self.log_factorial_map = log_factorial_map
        self.medium_pep_length = medium_pep_length
        #self.deleterious_pairs_record = deleterious_pairs_record
        
        
    # == Function Summary
    #   One-iteration algorithm to find one significant residue/position pair
    #   build the peptide motif from the iterative search
    #   Parameters: 
    #   -- foreground_peps_set, peptide dataset as foreground 
    #   -- self.fg_data_set, peptide dataset as background
    #   -- pseudo_count_matrix, the probability matrix generated from the beginning background set
    #   -- num_rp  - number of rp has been identified for the current motif
    #   -- sig_pairs - the pairs has been identified for the current motif
    #   Return:
    #   function return a tuple of values including Flag, motif_aa, motif_pos, logp_sig_pair, fg_vs_branch, bg_vs_branch
    #   example: (1, Y, 8, 1e-10, (25,100), (50,200))
    #   
    def find_one_pair(self, deleterious_pairs_record): 
        
        if len(self.bg_data_set) == 0:
            Flag = 2
            print("Background does not contain this motif, exit from search")
            return (Flag, '137', 137, 0)   
        
        Flag =0 # initialize flag
        
        #----------------------------------Step 1, build matrices --------------------------------------
        # Build the foreground_count_matrix in both forward and backward counting using dictionary type
        fw_start_time = time.time() 
        fg_matrix = bm.Matrix(self.MAX_LENGTH, self.fg_data_set) 
             
        fg_NT_char = fg_matrix.charMatrixAlignedNT()
        NT_fg_count_matrix = fg_matrix.aaCountMatrix(fg_NT_char)
        
        fg_CT_char = fg_matrix.charMatrixAlignedCT()
        CT_fg_count_matrix = fg_matrix.aaCountMatrix(fg_CT_char)
        
        fw_end_time = time.time()
        fw_matrix_runtime = (fw_end_time - fw_start_time)/60
        #print "fg_matrix count time = ", fw_matrix_runtime, "min"
        
        # Build the background_probablity_matrix as a database using dictionary type
        bg_start_time = time.time()
        bg_matrix = bm.bgMatrixWithPseudoCount(self.MAX_LENGTH, self.bg_data_set, 
                                               self.NT_pseudo_count_matrix, self.CT_pseudo_count_matrix)
        NT_bg_freq_matrix = bg_matrix.NTmatrix()
        CT_bg_freq_matrix = bg_matrix.CTmatrix()
        
        bg_end_time = time.time()
        bg_matrix_runtime = (bg_end_time - bg_start_time)/60
        #print "bg_matrix probability matrix time = ", bg_matrix_runtime, "min"
                                                          
        # Build the binomial_probablity_matrix using dictionary type for both forward and backward foreground alignment
        logbp_start_time = time.time()
        
        NT_bpl_matrix = fg_matrix.logbpMatrix(NT_fg_count_matrix, NT_bg_freq_matrix, self.log_factorial_map)
        CT_bpl_matrix = fg_matrix.logbpMatrix(CT_fg_count_matrix, CT_bg_freq_matrix, self.log_factorial_map)
 
  
        logbp_end_time = time.time()
        logbp_runtime = (logbp_end_time - logbp_start_time)/60
        #print "fw_logbp matrix time = ", logbp_runtime, "min"
             
        NT_sig_pair_list = SigPair.sigPairList(self, NT_bpl_matrix)
        # Significant pair position as positive int for CT now
        CT_sig_pair_list_P = SigPair.sigPairList(self, CT_bpl_matrix)
        CT_sig_pair_list = [dict() for i in range(len(CT_sig_pair_list_P))]
        for i in range(len(CT_sig_pair_list_P)):
            CT_sig_pair_list[i]['pair'] = (CT_sig_pair_list_P[i]['pair'][0], CT_sig_pair_list_P[i]['pair'][1] - self.MAX_LENGTH) 
            CT_sig_pair_list[i]['log_prob'] = CT_sig_pair_list_P[i]['log_prob']          
        # for debugging, check foreground matrix and pair list
        #print "fg_peps_set = ", self.fg_data_set
        #print "NT_fg_count_matrix = ", NT_fg_count_matrix
        #print "NT_sig_pair_list =", NT_sig_pair_list
        #print "CT_fg_count_matrix = ", CT_fg_count_matrix
        #print "CT_sig_pair_list =", CT_sig_pair_list
        
        # Check the significant pair found in the iterative search             
        if len(NT_sig_pair_list) == 0 and len(CT_sig_pair_list) == 0:
            # No more significant residue/position pairs are discovered
            Flag = 1
            return (Flag, '9', 137, 0)
        # check if either list is empty
        elif len(NT_sig_pair_list) == 0:
            # Significant pair search only need to conduct in one list
            rp_aa, rp_pos, rp_logbp = SigPair.findSigPairInOneList(self, CT_sig_pair_list, CT_fg_count_matrix, 
                                                                   deleterious_pairs_record)
            if rp_aa == 'empty':
                Flag = 1
                return (Flag, '9', 137, 0)
        elif len(CT_sig_pair_list) == 0 or CT_sig_pair_list == None:
            # Significant pair search only need to conduct in one list
            rp_aa, rp_pos, rp_logbp = SigPair.findSigPairInOneList(self, NT_sig_pair_list, NT_fg_count_matrix,
                                                                   deleterious_pairs_record)
            if rp_aa == 'empty':
                Flag = 1
                return (Flag, '9', 137, 0)
        else: 
            rp_aa, rp_pos, rp_logbp = SigPair.findSigPairInBothLists(self, NT_sig_pair_list, 
                                      CT_sig_pair_list, NT_fg_count_matrix, CT_fg_count_matrix,
                                      deleterious_pairs_record) 
            if rp_aa == 'empty':
                Flag = 1
                return (Flag, '9', 137, 0)  
            
        self.mypair = (rp_aa, rp_pos)  # Get the significant pair with the best bp value and position combination
        # print out the significant bp value for debug file   
        print("log form of binomial probability of significant rp pair = "), rp_logbp
        
        # Check if using Z-score or not
        if self.Z_flag == 0:
            print "No Z-score used"
            return (Flag, rp_aa, rp_pos, rp_logbp)
        
        else:
                
            Zscore = checkZScore.CheckZScore(self.MAX_LENGTH, self.fg_data_set, self.iPairSet, NT_bg_freq_matrix, 
                                  CT_bg_freq_matrix, rp_logbp, self.NT_pseudo_count_matrix, 
                                  self.NT_pseudo_count_matrix, self.log_factorial_map)   
      
            z_score = Zscore.calculate_ZScore()
            print "this rp z_score = ", z_score
            if z_score > self.Z_flag:
                self.mypair = (rp_aa, rp_pos)
                # print out the significant bp value for debug file   
                print("log form of binomial probability of significant rp pair = "), rp_logbp
            else:
                print("This pair doesn't pass z_score test, return this motif")
                Flag = 1
                return (Flag, '9', 137, 0)
                    
        return (Flag, rp_aa, rp_pos, rp_logbp)
    
    def getSigPair(self):
        '''
        return the identified pair in a format of tuple
        '''
        return self.mypair
    
    # Traverse the binomial_log_matrix to determine the best hits based on the probability less 
    # than the threshold. Data structure used is a dictionary with a tuple (aa, pos) as key and 
    # the binomial probability as value
    def sigPairList(self, bplog_matrix):
        
        sig_pair_list = []
        for aa in AA_LIST:
            for pos in range(self.MAX_LENGTH):
                sig_pair_dict = {}
                # Check if the binomial probability == 0, this caused by either p == 0, which is background have 0 r/p pair
                # or p == 1, which is the case after pruning
                if (bplog_matrix[aa][pos] < self.BINOMIAL_PROB_THRESHOLD and bplog_matrix[aa][pos] != 0):
                    sig_pair = (aa, pos)
                    if sig_pair not in self.iPairSet:
                        sig_pair_dict['pair'] = sig_pair
                        sig_pair_dict['log_prob'] = bplog_matrix[aa][pos]
                        sig_pair_list.append(sig_pair_dict)
                        
        return sig_pair_list
    
    # Find significant pair in one list if another is empty
    def findSigPairInOneList(self, sig_pair_list, fg_count_matrix, deleterious_pairs_record):
    
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
                #Check if cross the mid-point
                if (sorted_sig_pair_list[0]['pair'][0] > 0.5*self.medium_pep_length+1 and
                    sorted_sig_pair_list[1]['pair'][0] < 0.5*self.medium_pep_length+1 and
                    utility.diffByPercent(sorted_sig_pair_list[0]['pair'][1], 
                                              sorted_sig_pair_list[1]['pair'][1]) < 0.05):
                        rp_aa = sorted_sig_pair_list[1]['pair'][0]
                        rp_pos = sorted_sig_pair_list[1]['pair'][1]
                        rp_logbp = sorted_sig_pair_list[1]['log_prob']
                else:
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
        
        
    # Find significant pair in one list if another is empty
    def findSigPairInBothLists(self, NT_sig_pair_list, CT_sig_pair_list, fw_fg_count_matrix, bw_fg_count_matrix, 
                               deleterious_pairs_record):
        
        #rp_aa, rp_pos, rp_logbp = ('debug', 'debug', 'debug')
        
        fw_rp_aa, fw_rp_pos, fw_rp_logbp = SigPair.findSigPairInOneList(self, NT_sig_pair_list, fw_fg_count_matrix,
                                                                        deleterious_pairs_record)
        bw_rp_aa, bw_rp_pos, bw_rp_logbp = SigPair.findSigPairInOneList(self, CT_sig_pair_list, bw_fg_count_matrix,
                                                                        deleterious_pairs_record)
        if fw_rp_aa == 'empty' and bw_rp_aa != 'empty':
            rp_aa = bw_rp_aa
            rp_pos = bw_rp_pos
            rp_logbp = bw_rp_logbp
            return rp_aa, rp_pos, rp_logbp
        elif fw_rp_aa != 'empty' and bw_rp_aa == 'empty':
            rp_aa = fw_rp_aa
            rp_pos = fw_rp_pos
            rp_logbp = fw_rp_logbp
            return rp_aa, rp_pos, rp_logbp
        elif fw_rp_aa == 'empty' and bw_rp_aa == 'empty':
            return 'empty', 'empty', 'empty'
        else:
            fw_sig_pair = (fw_rp_aa, fw_rp_pos)
            bw_sig_pair = (bw_rp_aa, bw_rp_pos)
        
            seqMatchedFW = SigPair.seqSetMatchedPair(self, fw_sig_pair, self.fg_data_set)
            seqMatchedBW = SigPair.seqSetMatchedPair(self, bw_sig_pair, self.fg_data_set)
        
            #print "seqMatchedFW = ", seqMatchedFW
            print "len(seqMatchedFW)= ", len(seqMatchedFW)
            #print "seqMatchedBW = ", seqMatchedBW
            print "len(seqMatchedBW)= ", len(seqMatchedBW)
            print "FW pair, log_bp = ", fw_sig_pair, fw_rp_logbp
            print "BW pair, log_bp = ", bw_sig_pair, bw_rp_logbp
            
            #check if they are same
            if SigPair.isSameRP(self, fw_sig_pair, bw_sig_pair, seqMatchedFW, seqMatchedBW):
            
                if len(seqMatchedFW) >= len(seqMatchedBW):
                    if (fw_rp_pos > (self.medium_pep_length/2 + 1) and 
                        abs(bw_rp_pos) < (self.medium_pep_length/2 + 1) and
                        utility.diffByPercent(fw_rp_logbp, bw_rp_logbp) < 0.25):
                            rp_aa = bw_rp_aa
                            rp_pos = bw_rp_pos
                            rp_logbp = bw_rp_logbp
                    else:
                        rp_aa = fw_rp_aa
                        rp_pos = fw_rp_pos
                        rp_logbp = fw_rp_logbp
                else:
                    if (abs(bw_rp_pos) > (self.medium_pep_length/2 + 1) and 
                        fw_rp_pos < (self.medium_pep_length/2 + 1) and
                        utility.diffByPercent(fw_rp_logbp, bw_rp_logbp) < 0.25):
                            rp_aa = fw_rp_aa
                            rp_pos = fw_rp_pos
                            rp_logbp = fw_rp_logbp
                    else:
                        rp_aa = bw_rp_aa
                        rp_pos = bw_rp_pos
                        rp_logbp = bw_rp_logbp
            else:
                if fw_rp_logbp <= bw_rp_logbp:
                    if (fw_rp_pos > (self.medium_pep_length/2 + 1) and 
                        abs(bw_rp_pos) < (self.medium_pep_length/2 + 1) and
                         utility.diffByPercent(fw_rp_logbp, bw_rp_logbp) < 0.1):
                        rp_aa = bw_rp_aa
                        rp_pos = bw_rp_pos
                        rp_logbp = bw_rp_logbp
                        print "pick the backward position since diff = ", utility.diffByPercent(fw_rp_logbp, bw_rp_logbp)
                    else:
                        rp_aa = fw_rp_aa
                        rp_pos = fw_rp_pos
                        rp_logbp = fw_rp_logbp
                else:
                    if (fw_rp_pos < (self.medium_pep_length/2 + 1) and 
                        abs(bw_rp_pos) > (self.medium_pep_length/2 + 1) and
                        utility.diffByPercent(fw_rp_logbp, bw_rp_logbp) < 0.1):
                            rp_aa = fw_rp_aa
                            rp_pos = fw_rp_pos
                            rp_logbp = fw_rp_logbp
                            print "pick forward position since diff = ", utility.diffByPercent(fw_rp_logbp, bw_rp_logbp)
                    else:
                        rp_aa = bw_rp_aa
                        rp_pos = bw_rp_pos
                        rp_logbp = bw_rp_logbp
                     
            return rp_aa, rp_pos, rp_logbp
    
    # This function finds if the backwards significant pair is the same as forward significant pair
    # rp_pair1, rp_pair2 = rp pair, example ('R', 0), ('R', -9)
    def isSameRP(self, rp_fw, rp_bw, fw_set_match, bw_set_match):
        return SigPair.isSameAA(self, rp_fw[0], rp_bw[0]) and SigPair.isOverlap(self, fw_set_match, bw_set_match)  
    
    # compare two amino acids and return true if they are equal        
    def isSameAA(self, aa1, aa2):
        print "aa1=", aa1
        print 'aa2=', aa2
        return aa1 == aa2
    
    # Check the percent of overlap of two pairs in order to determine if they are actually the same pair
    def isOverlap(self, fw_set_match, bw_set_match):
        intersection = fw_set_match & bw_set_match
        print "intersection = ", intersection
        return len(intersection)/float(len(fw_set_match)) > 0.5 or len(intersection)/float(len(bw_set_match)) > 0.5
    
    # Get the number matched from both forward and backward matrix
    def seqSetMatchedPair(self, rp, data_set):
        seqSet = set([])
        for seq in data_set:                
            if len(seq) <= rp[1] or len(seq) < -(rp[1]):
                continue
            else:
                if rp[0] == seq[rp[1]]:
                    seqSet.add(seq)
        return seqSet
