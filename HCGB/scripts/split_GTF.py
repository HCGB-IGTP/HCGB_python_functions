#!/usr/bin/env python3
#############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy              ##
## Copyright (C):2019-2021 Lauro Sumoy Lab, IGTP, Spain    ##
#############################################################
"""
split_GTF converts GTF file into many small GTF files.
Usage: split_GTF.py {OPTIONS} [.GTF file]

It splits GTF into given number of files. It takes into account no genes or transcript are broken.
It is also possible to split according to chromosome (one gtf/chromosome)
"""

import os
import re
import argparse;
import traceback
from termcolor import colored

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.files_functions as HCGB_files

############################################################
def create_names_GTF(GTFfile, name_file, chr_option, num_files, debug=False):
    """
    Creates names for subsets of GTF files generated.
    
    It retrieves list of sequences in the GTF, and generates name for each subset.
    """
    ##
    dict_files_generated = {}
    
    if debug:
        print()
        HCGB_aes.debug_message("********************** create_names_GTF ********************** ", "yellow")
        HCGB_aes.debug_message("Create name for GTF subsets", "yellow")
    
    ## Chromosome option
    if chr_option:
        
        ## get list of entries
        with open(GTFfile) as f:
            list2 = [row.split()[0] for row in f]
            
            ## get uniq
            list2  = set(list2)
            
            ## remove # characters
            list2 = [x for x in list2 if not x.startswith('#')]
            
        
        ## create files for list of entries
        for seq in list2:
            file_name = name_file + "-Chr_" + str(seq) + ".gtf"
            dict_files_generated["Chr_" + str(seq)] = file_name

    if debug:
        print()
        HCGB_aes.debug_message("dict_files_generated: " + str(dict_files_generated), "yellow")
        print()
        HCGB_aes.debug_message("--------------------- create_names_GTF --------------------- ", "yellow")
        print()
        
    ## return
    return dict_files_generated

############################################################
def split_GTF_call(GTFfile, num_files, name, chr_option, path_given=False, debug=False):
    """
    This functions checks if it has been done previously the split of GTF. 
    If done, returns dict with files names generated using create_names_GTF.
    If GTF split not done, it calls split_GTF function. 
    
    :param GTFfile: Absolute path to GTF file to split
    :param num_files: Number of files to create with equal number of lines
    :param name: Name to include in the files names generated. By default, GTF basename included.
    :param chr_option: TRUE/FALSE If chr_option provided, split GTF into chromosome, scaffolds or reference sequences.
    :param path_given: Path to save results. Default use absolute path of GTF file provided
    :param debug: TRUE/FALSE for debugging messages
    """

    ## debug messaging    
    if debug:
        print()
        HCGB_aes.debug_message("********************** split_GTF_call ********************** ")
        HCGB_aes.debug_message("Checking if GTF has been previously splitted", "yellow")
    
    ## get absolute path and name
    name_file = get_path_name(GTFfile, path_given, name, debug=debug)

    filename_stamp = path_given + '/.split_GTF_success'
    if os.path.isfile(filename_stamp):
        
        if debug:
            HCGB_aes.debug_message("Timestamp exists: .split_GTF_success ", "yellow")

        ## check file names generated and return names
        files_generated = create_names_GTF(GTFfile, name_file, chr_option, num_files, debug=debug)
        
        re_run=False
        for files_GTF in files_generated.values():
            if not HCGB_files.is_non_zero_file(files_GTF):
                re_run=True
                
                if debug:
                    print()
                    HCGB_aes.debug_message("It is required to re-run split GTF: File is zero or does not exists\n" + files_GTF, "yellow")
                
                break
        if not re_run:
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, 'split GTF'), 'yellow'))
            return (files_generated)

    ## call to split GTF
    files_generated = split_GTF(GTFfile, num_files, name, chr_option, path_given, debug)
    
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)
    
    if debug:
        HCGB_aes.debug_message("********************** split_GTF_call ********************** ")
    
    
    return (files_generated)    
    
############################################################
def get_path_name(GTFfile, path_given, name, debug=False):
    
    ## get absolute path
    if (path_given):
        ## save files in given dir
        path_given = os.path.abspath(path_given)
    else:
        path_given = os.path.dirname(GTFfile)

    ## get name
    if not name:
        name=HCGB_files.get_file_name(GTFfile)

    ## crete name and path
    name2 = os.path.join(path_given, name)
    
    ## debugging messages
    if debug:
        print()
        HCGB_aes.debug_message("********************************** get_path_name **********************************", "yellow")
        HCGB_aes.debug_message("GTFfile: " + GTFfile, "yellow")
        HCGB_aes.debug_message("path_given: " + path_given, "yellow")
        HCGB_aes.debug_message("name: " + name, "yellow")        
        HCGB_aes.debug_message("name2: " + name2, "yellow")
        HCGB_aes.debug_message("--------------------- get_path_name --------------------- ", "yellow")
    ## 
    return (name2)

############################################################
def split_GTF(GTFfile, num_files, name, chr_option, path_given=False, debug=False):
    """
    This functions splits given GTF into multiple files, either a given number of files or
    one for each chromosome.

    :param GTFfile: Absolute path to GTF file to split
    :param num_files: Number of files to create with equal number of lines
    :param name: Name to include in the files names generated. By default, GTF basename included.
    :param chr_option: TRUE/FALSE If chr_option provided, split GTF into chromosome, scaffolds or reference sequences.
    :param path_given: Path to save results. Default use absolute path of GTF file provided
    :param debug: TRUE/FALSE for debugging messages
    """
    ## init dict to store files generated
    dict_files_generated = {}

    if debug:
        print()
        HCGB_aes.debug_message("*************************** split_GTF ***************************", "yellow")
        HCGB_aes.debug_message("GTFfile: " + GTFfile, "yellow")
        HCGB_aes.debug_message("num_files: " + str(num_files), "yellow")        
        HCGB_aes.debug_message("chr_option: " + str(chr), "yellow")

    ## Check file is readable
    if not (HCGB_files.is_non_zero_file(GTFfile)):
        print("ERROR: File not readable. Please check path for file:\n" + GTFfile)
        exit()

    ## get absolute path and name
    name = get_path_name(GTFfile, path_given, name, debug=debug)

    if debug:
        HCGB_aes.debug_message("path_given: " + path_given, "yellow")
        HCGB_aes.debug_message("name: " + name, "yellow")
    
    print("")
        
    ## Get options
    if (chr_option):
        ## Prevalence of chr split if provided.
        print("+ Splitting file by sequence...")
        
        try:
            #read a file
            fileReader = open(GTFfile)
            
            ## skip comments at the beginning of files
            while True:
                line = fileReader.readline()
                if not line.startswith('#'):
                    break
                
            lineCount=0
            line = fileReader.readline()
            field=line.strip().split('\t')
            chrid=field[0]

            file_name = name + "-Chr_" + str(chrid) + ".gtf"
            dict_files_generated["Chr_" + str(chrid)] = file_name
             
            stop=False
            line2=""

            while line != '':#empty is EOF
                ## create new file
                if lineCount == 0:
                    #create a file in append mode
                    field=line.strip().split('\t')
                    chrid=field[0]

                    file_name = name + "-Chr_" + str(chrid) + ".gtf"
                    dict_files_generated["Chr_" + str(chrid)] = file_name
                    fileWriter = open(file_name,"a") ## append as it might be some repetiteve elements at the end                

                #write a line
                fileWriter.write(line)
                
                ## Debug messages
                if debug:
                    HCGB_aes.debug_message("Chr: " + str(chrid), "red")
                    HCGB_aes.debug_message("line: " + line, "red")
                
                ## stop when Chr changes
                while True:
                    field=line.strip().split('\t')
                    chrid=field[0]

                    ## read new line
                    line2 = fileReader.readline()
                    field2=line2.strip().split('\t')
                    chrid2=field2[0]

                    ## debug messages    
                    if debug:
                        HCGB_aes.debug_message("line: " + line, "red")
                        HCGB_aes.debug_message("line2: " + line2, "red")
                        HCGB_aes.debug_message("geneid: " + chrid, "yellow")
                        HCGB_aes.debug_message("chrid2: " + chrid2, "yellow")

                    ##
                    if (chrid == chrid2):
                        fileWriter.write(line2)
                        line=line2
                    else:
                        ## init
                        lineCount = 0
                        fileWriter.close()
                        line=line2
                        break

            ## Close GTF File        
            fileWriter.close()
        
        except Exception as e:
            #print the exception if any
            print(e.__traceback__)
            traceback.print_exc()
        finally:
            #close the file reader
            fileReader.close()

    else:
        ## Split into several files as provided.

        # max lines you want to write in a single file
        totalLines_file = sum(1 for line in open(GTFfile))
        fileLineCount = int(totalLines_file/num_files)
        lineCount = 0
        fileCount = 1    

        ## Debug messages
        if debug:
            HCGB_aes.debug_message("totalLines_file: " + str(totalLines_file), "yellow")
            HCGB_aes.debug_message("fileLineCount: " + str(fileLineCount), "yellow")

        try:
            #read a file
            fileReader = open(GTFfile)
            line = fileReader.readline()
            file_name = name + "-" + str(fileCount) + ".gtf"
            dict_files_generated[name + "-" + str(fileCount)] = file_name

            stop=False
            line2=""

            while line != '':#empty is EOF
                ## create new file
                if lineCount == 0:
                    #create a file in append mode
                    file_name = name + "-" + str(fileCount) + ".gtf" 
                    dict_files_generated[name + "-" + str(fileCount)] = file_name
                    fileWriter = open(file_name,"w")
                    #increment file count, use it for new file name
                    fileCount += 1

                #write a line
                fileWriter.write(line)
                
                ## sum to lines
                lineCount += 1

                ## Debug messages
                if debug:
                    HCGB_aes.debug_message("lineCount: " + str(lineCount), "red")
                    HCGB_aes.debug_message("line: " + line, "red")
                
                ## stop when max_lines_file achieved
                if lineCount == fileLineCount:
                    ## control we are not splitting genes
                    while True:
                        field=line.strip().split('\t')
                        geneid=re.findall(r'gene_id \"([\w\.]+)\"',field[8])

                        ## read new line
                        line2 = fileReader.readline()
                        lineCount2 = lineCount + 1
                        field2=line2.strip().split('\t')
                        geneid2=re.findall(r'gene_id \"([\w\.]+)\"',field2[8])

                        ## debug messages    
                        if debug:
                            HCGB_aes.debug_message("lineCount: " + str(lineCount), "red")
                            HCGB_aes.debug_message("line: " + line, "red")
                            HCGB_aes.debug_message("lineCount2: " + str(lineCount2), "red")
                            HCGB_aes.debug_message("line2: " + line2, "red")
                            HCGB_aes.debug_message("geneid: " + geneid[0], "yellow")
                            HCGB_aes.debug_message("geneid2: " + geneid2[0], "yellow")
                        
                        ##
                        if (geneid[0] == geneid2[0]):
                            fileWriter.write(line2)
                        else:
                            stop=True
                            break

                    ## init
                    lineCount = 0
                    fileWriter.close()
                
                if not stop:
                    #read a line
                    line = fileReader.readline()
                    if line == '':#empty is EOF
                        fileWriter.close()
                else:
                    line=line2
                    stop=False
            
            ## Close GTF File        
            fileWriter.close()
        
        except Exception as e:
            #print the exception if any
            print(e.__traceback__)
            traceback.print_exc()
        finally:
            #close the file reader
            fileReader.close()

    ##
    if debug:
        print()
        HCGB_aes.debug_message("dict_files_generated: " + str(dict_files_generated), "yellow")
        print()
        HCGB_aes.debug_message("*************************** split_GTF ***************************", "yellow")
        
    
    return(dict_files_generated)

############################################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Split GTF file into multiple files
    Note:
    It takes into account no genes or transcript are broken
    It is also possible to split according to chromosome (one gtf/chromosome)
    
    ''');
    
    parser.add_argument('--input', '-i', help='Input GTF file', required=True);
    
    parser.add_argument('--num_files','-n', type=int,
                        help='Split GTF file into as many subfiles.', default=2);
    
    parser.add_argument('--name',
                        help='Name to add for each file generated. Default: use filename provided.', default="");

    parser.add_argument('--path',
                        help='Path to save for each file generated. Default: use path from filename provided.', default="");

    parser.add_argument('--split_chromosome','-c',action="store_true",
                        help='Split GTF file for each chromosome or reference sequence available.');

    args=parser.parse_args();
    
    ## lets split the big file provided
    files_generated = split_GTF(os.path.abspath(args.input), num_files=args.num_files, name=args.name, 
              chr_option=args.split_chromosome, path_given=args.path, debug=True)
    
    print(files_generated)
    

############################################################
if __name__== "__main__":
    main()