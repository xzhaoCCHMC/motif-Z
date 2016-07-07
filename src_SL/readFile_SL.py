# This module handles file reading
# -- csv file from gfy database searching 
# -- csv/txt file of extracted raw peptide in gfy format 

import csv
import os
import re


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
    if os.path.exists(csv_file):
        output_file = csv_file.split(".")[0]
        fw = open(output_file+".dat", 'w')
        count = 0
        pep = set([])
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
                pep.add(peptide)
                count=count+1
                
        print("Total sequence in " + csv_file + ": " + str(count))
        print("Total unique peptides in " + csv_file + ": " + str(len(pep)))
        fw.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep
 

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
    if os.path.exists(txt_file):
        count = 0
        pep = set([])
        with open(txt_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                row = row.rstrip('\n')
                #print(row)
                pep.add(row)
                count=count+1
                
        print("Total sequence in " + txt_file + ": " + str(count))
        print("Total unique peptides in " + txt_file + ": " + str(len(pep)))
        f.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep
   

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
                if row != []:
                    peptide = row[0]
                    fw.write(peptide + '\n')
                    pep.add(peptide)
                    count=count+1
                
        print("Total sequence in " + csv_file + ": " + str(count))
        print("Total unique peptides in " + csv_file + ": " + str(len(pep)))
        fw.close()
    else:
        print("The file is not exist, please check the path and input again")
    return pep



  




