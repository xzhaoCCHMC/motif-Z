# This module handles file reading
# -- csv file from gfy database searching 
# -- csv/txt file of extracted raw peptide in gfy format 

import csv
import os
import re
import utility


# == Function Summary
#   General function to take input file name and output the filename with path to the main function
#   Parameters:
#   -- prompt, a string to give the correct file path
#   Return:
#   --The file path input by user
def readin_filenmame(prompt):
    while True:
        motif_file = raw_input(prompt)
        if (os.path.exists(motif_file)):
            break
        else:
            print ("The file"  + motif_file + " doesn't exist.")
    return motif_file

def read_textInput(textInput): 
    peps_list= textInput.split('\n')
    peps_set = set([])
    data_summary = {}
    count=0
    for pep in peps_list:
        pep = pep.rstrip()
        count += 1
        peps_set.add(pep)
        length = len(pep)
        if not data_summary.has_key(length):
            data_summary[length] = 1
        else:
            data_summary[length] += 1
                  
    print("Total sequence in text input: " + str(count))
    print("Total unique peptides in text input" + str(len(peps_set))) 
    
    return peps_set, data_summary         
                
             
# == Function Summary
#   Read in .csv or .txt file 
#   the peptides data from LC-MS/MS run and SEQUEST search
#   function process input data to extract only peptide sequence between dots. 
#   Parameters:
#   -- csv_file, file name which contains the SEQUEST results
#   -- opt_dataSep, the separators/pattern in the input data file: "\." for raw gfy file of binding peptides; 
#   -- most parameters have default value for binding peptide sequence
#   Return:
#   -- Set which has the all unique peptides 
#   -- A text file in which all peptide sequences written in the one-peptide-per-line style
def read_gfyPep(csv_file, opt_sep = '\.'):
    pep = set([])
    data_summary = {}
    if os.path.exists(csv_file):
        output_file = csv_file.split(".")[0]
        fw = open(output_file+".dat", 'w')
        count = 0
        with open(csv_file, 'rU') as f:
            next(f)
            reader = csv.reader(f)
            
            for row in reader:
                #Check reading
                #print(row)
                if opt_sep == '\.':
                    if len(re.split(opt_sep, row[0])) == 2:
                        peptide = re.split(opt_sep, row[0])[0]
                    else:
                        peptide = re.split(opt_sep, row[0])[1]
                else:
                    print('check the format of peptides, not separated by dot')
                fw.write(peptide + '\n')
                length = len(peptide)
                if not data_summary.has_key(length):
                    data_summary[length] = 1
                else:
                    data_summary[length] += 1
                pep.add(peptide)
                count=count+1
                
        print("Total sequence in " + csv_file + ": " + str(count))
        print("Total unique peptides in " + csv_file + ": " + str(len(pep)))
        fw.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep, data_summary
 

# == Function Summary
#   Read in .csv or .txt file 
#   the peptides data from LC-MS/MS run and SEQUEST search
#   function process input data to extract only peptide sequence between dots. 
#   Parameters:
#   -- csv_file, file name which contains the SEQUEST results
#   -- opt_dataSep, the separators/pattern in the input data file: "\." for raw gfy file of binding peptides; 
#   -- most parameters have default value for binding peptide sequence
#   Return:
#   -- Set which has the all unique peptides 
#   -- A text file in which all peptide sequences written in the one-peptide-per-line style
def read_txtPep(txt_file):
    data_summary = {}
    if os.path.exists(txt_file):
        count = 0
        pep_set = set([])
        with open(txt_file, 'r') as f:
#            lines = f.readlines()
#            for row in lines:
#                row = row.rstrip('\n')
#                print(row)
#                pep.add(row)
#                count=count+1
            line = f.readline()
            while line:
                line = line.rstrip()
                count += 1
                pep_set.add(line)
                length = len(line)
                if not data_summary.has_key(length):
                    data_summary[length] = 1
                else:
                    data_summary[length] += 1
                line = f.readline()      
        print("Total sequence in " + txt_file + ": " + str(count))
        print("Total unique peptides in " + txt_file + ": " + str(len(pep_set)))
        f.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep_set, data_summary


# == Function Summary
#   Read in .csv file without Gfy format 
#   the peptides data from LC-MS/MS run and SEQUEST search
#   function process input data to extract only peptide sequence between dots. 
#   Parameters:
#   -- csv_file, file name which contains the SEQUEST results
#   -- opt_dataSep, the separators/pattern in the input data file: "\." for raw gfy file of binding peptides; 
#   -- most parameters have default value for binding peptide sequence
#   Return:
#   -- Set which has the all unique peptides 
#   -- A text file in which all peptide sequences written in the one-peptide-per-line style
def read_csvPep(csv_file):
    pep = set([])
    data_summary = {}
    if os.path.exists(csv_file):
        count = 0
        with open(csv_file, 'rU') as f:
            next(f)
            reader = csv.reader(f)
            
            for row in reader:
                #Check reading
                #print(row)
                if row != []:
                    peptide = row[0]
                    length = len(peptide)
                    if not data_summary.has_key(length):
                        data_summary[length] = 1
                    else:
                        data_summary[length] += 1
                    pep.add(peptide)
                    count=count+1
                
        print("Total sequence in " + csv_file + ": " + str(count))
        print("Total unique peptides in " + csv_file + ": " + str(len(pep)))
        f.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep, data_summary

# == Function Summary
#   Read in .txt result files
#   Return:
#   -- Set which has the all unique peptides 
#   -- A text file in which all peptide sequences written in the one-peptide-per-line style
def read_resultfile(result_file, MOTIF_LENGTH, MINIMUM_LENGTH):
    motifs = set([])
    if os.path.exists(result_file):
        count = 0
        with open(result_file, 'rU') as f:
            for i in range(5):
                next(f)
            for line in f:
#            while True:
#                line = f.readline()
                if not line: break
                motif = line.split('\t')[0]
                motif = motif.rstrip()
                print(motif)
                rexp_motif = utility.format_to_rexp(motif, MOTIF_LENGTH, MINIMUM_LENGTH)
                motifs.add(rexp_motif)
                count=count+1
            
        print("Total motifs in " + result_file + ": " + str(count))
        f.close()
    else:
        print("The file is not exist, please check the path and input again")
    return motifs




  




