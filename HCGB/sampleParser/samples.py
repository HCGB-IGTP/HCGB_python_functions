#!/usr/bin/ python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Prepares samples for further analysis.
"""

import os
import re
import pandas as pd
from termcolor import colored

def select_samples (list_samples, samples_prefix, pair=True, exclude=False, Debug=False, lane=False, include_all=False):
    """
    Select samples
    
    Given a sample prefix (any or a given list), this function retrieves
    sample files from a list given. If exclude option provided, excludes 
    the files retrieved from the total.
    
    :param list_samples: List of absolute path for fastq files
    :param samples_prefix: List of prefix to search 
    :param pair: True/false for paired-end files
    :param exclude: True/false for exclude found files from the total
    :param Debug: True/false for debugging messages
    :param lane: Include lane tag within name id
    
    :type list_samples: list
    :type samples_prefix: list
    :type pair: bool
    :type exclude: bool
    :type Debug: bool
    :type lane: bool
    
    :returns: Dataframe
    """
    #Get all files in the folder "path_to_samples"
    sample_list = pd.DataFrame(columns=('sample', 'file'))
    
    for names in samples_prefix:
        for path_fastq in list_samples:    
            fastq = os.path.basename(path_fastq)
            samplename_search = re.search(r"(%s)\_{0,1}(R1|1|R2|2){0,1}(.*){0,1}\.f.*q.*" % names, fastq)
            enter = ""
            if samplename_search:
                if (exclude): ## exclude==True
                    enter = False
                else: ## exclude==True
                    enter = True
            else:
                if (exclude): ## exclude==True
                    enter = True
                else: ## exclude==True
                    enter = False
                    
            if enter:
                if fastq.endswith('.gz') or fastq.endswith('fastq') or fastq.endswith('fq'):
                    sample_list.loc[len(sample_list)] = (names, path_fastq) 
                else:
                    ## debug message
                    if (Debug):
                        print (colored("**DEBUG: sampleParser.select_samples **", 'yellow'))
                        print (colored("** ERROR: %s is a file that is neither in fastq.gz or .fastq format, so it is not included" %path_fastq, 'yellow'))
                            
    ## discard duplicates if any
    non_duplicate_names = sample_list['sample'].to_list() #
    non_duplicate_names = list(set(non_duplicate_names))
    
    ## it might be a bug in exclude list.
    ## if sample X1 is provided to be excluded, we might be also excluding
    ## sample X12, sample X13, etc.
    ## TODO: check this

    ## debugging messages
    if Debug:
        print (colored("** DEBUG: select_samples",'yellow'))
        print ("non_duplicate_names:")
        print (non_duplicate_names)
    
    ## check they match with given input
    if (exclude): ## exclude==True
        if bool(set(samples_prefix).intersection(non_duplicate_names)):
            print(colored("** ERROR: Some non desired samples are included", 'red'))
    else: ## exclude==True
        non_duplicate_names = set(samples_prefix).intersection(non_duplicate_names)

    ## get fields
    
    tmp = sample_list[ sample_list['sample'].isin(non_duplicate_names) ]
    non_duplicate_samples = tmp['file'].to_list()
    
    ## debugging messages
    if Debug:
        print (colored("** DEBUG: select_samples",'yellow'))
        print ("non_duplicate_names:")
        print (non_duplicate_names)
        print ("samples_prefix")
        print (samples_prefix)
        print ("non_duplicate_samples")
        print (non_duplicate_samples)
        print ("tmp dataframe")
        #functions.print_all_pandaDF(tmp)
        print(tmp)
                
    ## get info
    name_frame_samples = get_fields(non_duplicate_samples, pair, Debug, include_all)    
    number_files = name_frame_samples.index.size
    total_samples = set(name_frame_samples['name'].to_list())
    
    ##
    if (lane):
        ## include lane tag within name
        name_frame_samples['name'] = name_frame_samples['name'] + '_' + name_frame_samples['lane']
        name_frame_samples['new_name'] = name_frame_samples['name']
            
    ## debugging messages
    if Debug:
        print (colored("** DEBUG: select_samples",'yellow'))
        print ("name_frame_samples:")
        print (name_frame_samples)
        print ("number_files:")
        print (number_files)
        print ("total_samples:")
        print (total_samples)
    
    ### get some stats
    if (number_files == 0):
        print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
        exit()
    print (colored("\t" + str(number_files) + " files selected...", 'yellow'))
    print (colored("\t" + str(len(total_samples)) + " samples selected...", 'yellow'))
    if (pair):
        print (colored("\tPaired-end mode selected...", 'yellow'))
    else:
        print (colored("\tSingle end mode selected...", 'yellow'))
    
    ## return info
    return (name_frame_samples)

###############
