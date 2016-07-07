# This module handles parsing command line parameters
# defines what arguments it requires, and argparse will figure out how to parse those out of sys.argv. 
# The argparse module also automatically generates help and usage messages and issues errors 
# when users give the program invalid arguments.
# -- take command line arguments and store them in the list args 
# -- print out argument for checking

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--ft", help="Input peptide file type, .csv with or without gfy format or .txt", dest="fileType",
                     default="csv", type=str, choices = ["csv", "csvGfy", "txt"])
parser.add_argument("--inputf", help="Input peptide file with correct path",metavar="FILE",
                     required=True, dest="inputFile")
parser.add_argument("--length", help="number of amino acids in the peptides",
                     default=9, type=int, choices=[7,8,9,10,11,12,23])
parser.add_argument("--bp", help="threshold of binomial probability of a residue/position pair",
                     default=1000, type=float)
parser.add_argument("--zscore", help="z-score threshold to accept the significant rp pair",
                     default='4', type=str)
parser.add_argument("--thresholdFG", help="the smallest number of peptides in foreground to be considered including the motif",
                     default=20, type=int)
parser.add_argument("--outputf", help="Output file name for result and debug without suffix", default='results/mySearch')
parser.add_argument("--verbose", help="Verbose output and stdio errors",
                     action="store_true")

#Get all the argument and store them in a global structure
args = parser.parse_args()

# /Users/zhaoxh/Documents/workspace/Motif-X_v1/readFile_test.csv