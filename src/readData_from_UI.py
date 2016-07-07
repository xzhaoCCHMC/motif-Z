# This module handles file reading and input text from UI 


import readFile as rf


# == Function Summary
#   General function to take input file name and output the filename with path to the main function
#   Parameters:
#   -- file_type: input from UI by user
#   -- pep_file, file path saved in model
#   Return:
#   --peptides in set and data_summary
def readFile_from_UI(file_type, pep_file):
    # Build foreground dataset from input peptide file
    if file_type == "csv":
        peps, data_summary = rf.read_csvPep(pep_file)
    if file_type == "txt":
        peps, data_summary = rf.read_txtPep(pep_file)
        
    return peps, data_summary

# == Function Summary
#   General function to take input file name and output the filename with path to the main function
#   Parameters:
#   -- prompt, a string to give the correct file path
#   Return:
#   --The file path input by user
def readInput_from_UI(pep_input):
    # Build foreground dataset from input peptide file
    peps, data_summary = rf.read_textInput(pep_input)
        
    return peps, data_summary
