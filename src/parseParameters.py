# This module handles parsing command line parameters
# defines what arguments it requires, and argparse will figure out how to parse those out of sys.argv. 
# The argparse module also automatically generates help and usage messages and issues errors 
# when users give the program invalid arguments.
# -- take command line arguments and store them in the list args 
# -- print out argument for checking

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--ft", help="Input peptide file type, .csv with or without gfy format or .txt", dest="fileType",
                     default="csv", type=str, choices = ["csv", "csvGFY", "txt"])
parser.add_argument("--inputf", help="Input peptide file name with correct path",metavar="FILE",
                     required=True, dest="inputFile")
parser.add_argument("--maxLength", help="length of the longest peptides in the analysis",
                     default=14, type=int, choices=[7,8,9,10,11,12,14])
parser.add_argument("--minLength", help="length of the shortest peptides in the analysis",
                     default=7, type=int, choices=[7,8,9,10,11,12,14])
parser.add_argument("--zscore", help="if use Z-score in the searching, input [1-6] ",
                     default="4", type=str)
parser.add_argument("--bp", help="Optional feature: threshold of binomial probability of a residue/position pair, e.g. 1e-6",
                     default=1000, type=float)
parser.add_argument("--thresholdFG", help="the smallest number of peptides in foreground to be considered including the motif, default = 20",
                     default=20, type=int)
parser.add_argument("--outputf", help="Output file name for all results as prefix with correct file path", default='results/mySearch',
                    required=True, dest="outputFile")
parser.add_argument("--verbose", help="Verbose output and stdio errors",
                     action="store_true")

#Get all the argument and store them in a global structure
args = parser.parse_args()

# /Users/zhaoxh/Documents/workspace/Motif-X_v1/readFile_test.csv