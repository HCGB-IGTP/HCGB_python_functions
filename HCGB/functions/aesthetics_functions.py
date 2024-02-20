#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Shared functions used along ``BacterialTyper`` & ``XICRA`` pipeline.
With different purposes:
    - Print time stamps
  
    - Create system calls

    - Manage/list/arrange files/folders

    - Aesthetics

    - Manage fasta files

    - Other miscellaneous functions
"""
## useful imports
import sys
from termcolor import colored

############################################################################
########                     AESTHETICS                             ########
############################################################################
def pipeline_header(option, ver=""):
    """
    Prints a common header for the pipeline including name, author, copyright and year.        
    """
    
    print ("\n")
    print_sepLine("#", 70, False)
    
    if (option == 'XICRA'):
        print('#', '{: ^66}'.format("XICRA pipeline"), '#')
        print('#', '{: ^66}'.format("Jose F. Sanchez & Lauro Sumoy"), '#')
    elif (option == 'BacterialTyper'):
        print('#', '{: ^66}'.format("BacterialTyper pipeline"), '#')
        print('#', '{: ^66}'.format("Jose F. Sanchez, Cristina Prat & Lauro Sumoy"), '#')
    elif (option == 'RSP'):
        print('#', '{: ^66}'.format("RNAseq pipeline - RSP"), '#')
        print('#', '{: ^66}'.format("Jose F. Sanchez, Mireia Marin & Lauro Sumoy"), '#')
    
    print('#', '{: ^66}'.format("Copyright (C) 2019-2024 Lauro Sumoy Lab, IGTP, Spain"), '#')
    
    ## add version
    if ver:
        print('#', '{: ^66}'.format("Version: " + str(ver)), '#')
    
    print_sepLine("#", 70, False)

###############
def print_sepLine(char, num, color):
    string = char * num
    if (color):
        print (colored(string, color))
    else:
        print (string)

###############
def boxymcboxface(message):
    ## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
    ## give credit to them appropiately
    #print('-' * 79)
    print ('\n')
    print('|', '=' * 50, '|', sep='')
    print('|', '{: ^48}'.format(message), '|')
    print('|', '=' * 50, '|', sep='')
    print ('\n')
    #print('-' * 79)

###############
def progbar(curr, total, full_progbar):
    frac = (curr/total)
    filled_progbar = round(frac*full_progbar)
    print ('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')
    sys.stdout.flush()
    
###############
def green_message(message, color='light_green'):
    print (colored(message, color))   

def debug_message(message, color='yellow'):
    print (colored("*** DEBUG: %s" %(message), color))   

###############
def error_message(message, color='light_red'):
    print (colored("*** ERROR: %s" %(message), color))   

def warning_message(message, color='cyan'):
    print (colored("*** WARNING: %s" %(message), color))   

def raise_and_exit(message):
    error_message(message, color='light_red')
    error_message("Exit BacterialTyper", color='light_red')
    raise SystemExit()
    
def print_argparse_dict(dict_options):
    debug_message("Options provided: ")
    print("\t-")
    print('\n\t- '.join(f'{k}={v}' for k, v in vars(dict_options).items()))