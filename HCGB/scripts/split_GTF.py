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
import sys
import re
import pandas as pd
import HCGB.functions.aesthetics_functions as HCGB_aes
import argparse;
import traceback
import HCGB.functions.files_functions as HCGB_files

############################################################
def split_GTF(GTFfile, num_files, name, chr, path_given=False, debug=False):
    """
    This functions splits given GTF into multiple files, either a given number of files or
    one for each chromosome.

    :param GTFfile: Absolute path to GTF file to split
    :param num_files: Number of files to create with equal number of lines
    :param name: Name to include in the files names generated. By default, GTF basename included.
    :param chr: TRUE/FALSE If chr option provided, split GTF into chromosome, scaffolds or reference sequences.
    :param path_given: Path to save results. Default use absolute path of GTF file provided
    :param debug: TRUE/FALSE for debugging messages
    """

    dict_files_generated = {}

    if debug:
        HCGB_aes.debug_message("split_GTF function", "yellow")
        HCGB_aes.debug_message("GTFfile: " + GTFfile, "yellow")
        HCGB_aes.debug_message("num_files: " + str(num_files), "yellow")        
        HCGB_aes.debug_message("name provided: " + name, "yellow")
        HCGB_aes.debug_message("chr: " + str(chr), "yellow")

    ## Check file is readable
    if not (HCGB_files.is_non_zero_file(GTFfile)):
        print("ERROR: File not readable. Please check path for file:\n" + GTFfile)
        exit()

    ## get absolute path
    if (path_given):
        ## save files in given dir
        path_given = os.path.abspath(path_given)
    else:
        path_given = os.path.dirname(GTFfile)

    ## get name
    if not name:
        name=HCGB_files.get_file_name(GTFfile)
        print(name)
        
    name = os.path.join(path_given, name)

    if debug:
        HCGB_aes.debug_message("path_given: " + path_given, "yellow")
        HCGB_aes.debug_message("name: " + name, "yellow")
    
    print("")
        
    ## Get options
    if (chr):
        ## Prevalence of chr split if provided.
        print("+ Splitting file by sequence...")
        
        try:
            #read a file
            lineCount=0
            fileReader = open(GTFfile)
            line = fileReader.readline()
            field=line.strip().split('\t')
            chrid=field[0]

            file_name = name + "-Chr_" + str(chrid) + ".txt"
            dict_files_generated["Chr_" + str(chrid)] = file_name
             
            stop=False
            line2=""

            while line != '':#empty is EOF
                ## create new file
                if lineCount == 0:
                    #create a file in append mode
                    field=line.strip().split('\t')
                    chrid=field[0]

                    file_name = name + "-Chr_" + str(chrid) + ".txt"
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
            file_name = name + "-" + str(fileCount) + ".txt"
            dict_files_generated[name + "-" + str(fileCount)] = file_name

            stop=False
            line2=""

            while line != '':#empty is EOF
                ## create new file
                if lineCount == 0:
                    #create a file in append mode
                    file_name = name + "-" + str(fileCount) + ".txt" 
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
              chr=args.split_chromosome, path_given=args.path, debug=True)
    

############################################################
if __name__== "__main__":
    main()