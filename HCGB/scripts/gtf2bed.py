#!/usr/bin/env python3
#############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy              ##
## Copyright (C):2012 Davud Li Wei - https://weililab.org/ ##
##               2019-2021 Lauro Sumoy Lab, IGTP, Spain    ##
#############################################################
"""
gtf2bed.py converts GTF file to BED file.
Usage: gtf2bed.py {OPTIONS} [.GTF file]
History
    Nov.5th 2012:
        1. Allow conversion from general GTF files (instead of only Cufflinks supports).
        2. If multiple identical transcript_id exist, transcript_id will be appended a string like "_DUP#" to separate.

    Nov.9th 2021: Modified by JFSanchezherrero
        1. Reformat code to use as a function
        2. Add additional details
        3. Simplify

Copyrigth: http://alumni.cs.ucr.edu/~liw/scripts.html
BED Format details: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
"""

import sys
import re
import pandas as pd
import HCGB.functions.aesthetics_functions as HCGB_aes

global allids
allids = {}

global data_bed
data_bed = pd.DataFrame(columns=["chr", "start", "end", "transcript_id", "gene_id", "transcript_biotype",
                                 "fpkm", "strand", "start", "end", "RGB", 
                                 "estart", "strs", "strl"])

############################################################
def printbedline(estart, eend, field, nline, debug=False):
    ## debug messages    
    if debug:
        HCGB_aes.debug_message("GTF line: ", "yellow")
        print(field)
        
    estp=estart[0]-1
    eedp=eend[-1]
    
    ## ---------------------------------
    # use regular expression to get transcript_id, gene_id and expression level
    ## ---------------------------------
    geneid=re.findall(r'gene_id \"([\w\.]+)\"',field[8])
    transid=re.findall(r'transcript_id \"([\w\.]+)\"',field[8])
    fpkmval=re.findall(r'FPKM \"([\d\.]+)\"',field[8])
    biotype=re.findall(r'transcript_biotype \"([\w\.]+)\"',field[8])
    ## add others if required: e.g. transcript_biotype
    
    ## ---------------------------------
    ## Get Gene ID
    ## ---------------------------------
    if len(geneid)==0:
        print('Warning: no gene_id field in line ' + nline, file=sys.stderr)
        gene_id="none"
    else:
        gene_id=geneid[0]
    
    ## ---------------------------------
    ## Get transcript biotype
    ## ---------------------------------
    if len(biotype)==0:
        print('Warning: no transcript_biotype field in line ' + nline, file=sys.stderr)
        transcript_biotype="none"
    else:
        transcript_biotype=biotype[0]

    ## ---------------------------------
    ## FPKM field
    ## ---------------------------------
    fpkmint=100
    if len(fpkmval)>0: # Warning: no FPKM field
        fpkmval=fpkmval[0]
        fpkmint=round(float(fpkmval))

    ## ---------------------------------
    ## Get transcript ID
    ## ---------------------------------
    if len(transid)==0: # Warning: no transcript_id field
        transid='Trans_'+str(nline) ## set a new name for transcript if missing
        print('Warning: no transcript_id field in line ' + nline, file=sys.stderr)
        print('Warning: Generate new: ' + transid, file=sys.stderr)
    else:
        transid=transid[0]
    
    ## Previous transcript ID
    if transid in allids.keys():
        transid2=transid+'_DUP'+str(allids[transid])
        allids[transid]=allids[transid]+1
    ## New transcript ID
    else:
        allids[transid]=1
    
    ## ---------------------------------
    ## Get exon start and lengths
    ## ---------------------------------
    seglen=[eend[i]-estart[i]+1 for i in range(len(estart))]
    if (field[6]=="+"):
        segstart=[estart[i]-estart[0] for i in range(len(estart))]
    elif (field[6]=="-"):
        segstart=[estart[0]-estart[i] for i in range(len(estart))]

    strl=str(seglen[0])
    for i in range(1,len(seglen)):
        strl+=','+str(seglen[i])

    strs=str(segstart[0])
    for i in range(1,len(segstart)):
        strs+=','+str(segstart[i])

    ## debug messages    
    if debug:
        HCGB_aes.debug_message("estart: " + str(estart), "yellow")
        HCGB_aes.debug_message("seglen: " + str(seglen), "yellow")
        HCGB_aes.debug_message("segstart: " + str(segstart), "yellow")
        HCGB_aes.debug_message("strl: " + str(strl), "yellow")
        HCGB_aes.debug_message("strs: " + str(strs), "yellow")
    
    ## ---------------------------------
    ## Save data information
    ## ---------------------------------
    data_bed.loc[transid, "chr"] = field[0] 
    data_bed.loc[transid, "start"] = str(estp)
    data_bed.loc[transid, "end"] = str(eedp) 
    data_bed.loc[transid, "transcript_id"] = transid
    data_bed.loc[transid, "gene_id"] = gene_id
    data_bed.loc[transid, "transcript_biotype"] = transcript_biotype
    data_bed.loc[transid, "fpkm"] = str(fpkmint)
    data_bed.loc[transid, "strand"] = field[6]
    data_bed.loc[transid, "thickstart"] = str(estp)
    data_bed.loc[transid, "thickend"] = str(eedp) 
    data_bed.loc[transid, "RGB"] = "255,0,0"
    data_bed.loc[transid, "estart"] = len(estart)
    data_bed.loc[transid, "strl"] = strl
    data_bed.loc[transid, "strs"] = strs
    
    ## debug messages    
    if debug:        
        HCGB_aes.debug_message("data for transid: " + transid, "yellow")
        print(data_bed.loc[transid])    

############################################################
def parse_GTF(gtf_file, out_file, debug=False):
    ## Start the parsing of GTF

    ## Init variables 
    estart=[]
    eend=[]
    nline=0 # read lines one to one
    prevfield=[]
    prevtransid=''
    
    ## Loop through big GTF file
    for lines in open(gtf_file):
        field=lines.strip().split('\t')
        ## count lines
        nline=nline+1
        
        ## debug messages    
        if debug:        
            HCGB_aes.debug_message("Line: ", "yellow")
            print(lines)
            HCGB_aes.debug_message("nline: " + str(nline), "yellow")

        ## skip comment lines
        if field[0].startswith("#"): ## Comment line: skipping
            continue
                
        if len(field)<9:
            print('Error: the GTF should has at least 9 fields at line '+str(nline),file=sys.stderr)
            continue

        if field[1]!='Cufflinks':
            pass
            #print('Warning: the second field is expected to be \'Cufflinks\' at line '+str(nline),file=sys.stderr)
        
        if field[2]!='exon' and field[2] !='transcript':
            #print('Error: the third filed is expected to be \'exon\' or \'transcript\' at line '+str(nline),file=sys.stderr)
            continue
        
        ## ---------------------------------
        # use regular expression to get transcript_id
        transid=re.findall(r'transcript_id \"([\w\.]+)\"',field[8])
        if len(transid)>0:
            transid=transid[0]
        else:
            transid=''
            
        if field[2]=='transcript' or (prevtransid != '' and transid!='' and transid != prevtransid):
            #print('prev:'+prevtransid+', current:'+transid)
            # A new transcript record, write
            if len(estart)!=0:
                
                ## debug messages    
                if debug:        
                    HCGB_aes.debug_message("printbedline call")
                    HCGB_aes.debug_message("estart: " + str(estart))
                    HCGB_aes.debug_message("eend: " + str(eend))
                    HCGB_aes.debug_message("prevfield: " + str(prevfield))
                
                ## save record in bed format
                printbedline(estart, eend, prevfield, nline, debug)

            # Reset
            estart=[]
            eend=[]
        
        ## ---------------------------------
        prevfield=field
        prevtransid=transid
        if field[2]=='exon':
            try:  
                est=int(field[3])
                eed=int(field[4])
                estart+=[est]
                eend+=[eed]
            except ValueError:
                print('Error: non-number fields at line '+str(nline),file=sys.stderr)
            
    # the last record
    if len(estart)!=0:
        ## debug messages    
        if debug:        
            HCGB_aes.debug_message("printbedline call")
            HCGB_aes.debug_message("estart: " + str(estart))
            HCGB_aes.debug_message("eend: " + str(eend))
            HCGB_aes.debug_message("prevfield: " + str(prevfield))

        ## save record in bed format
        printbedline(estart, eend, field, nline, debug)

    ## ---------------------------------
    ## print information info file and or return
    ## ---------------------------------
    ## remove some unnecessary columns in data_bed

    ## ---------------------------------
    ## return when finished
    return(data_bed)

############################################################
def main():
    ## this code runs when call as a single script
    if len(sys.argv)<2:
        print('This script converts .GTF into .BED annotations.\n')
        print('Usage: gtf2bed [.GTF file] [OUT bed]\n')
        print('\nNote:')
        print('1\tOnly "exon" and "transcript" are recognized in the feature field (3rd field).')
        print('2\tIn the attribute list of .GTF file, the script tries to find "gene_id", "transcript_id" and "FPKM" attribute, and convert them as name and score field in .BED file.') 
        
        print('Author: Wei Li (li.david.wei AT gmail.com)')
        sys.exit()

    ## get output file
    out_file=""
    if len(sys.argv)>=3:
        out_file = sys.argv[2]
    else:
        out_file = "example.bed"
    
    ## parse 
    parse_GTF(sys.argv[1], out_file, False)

    ## final
    print(data_bed)


############################################################
if __name__== "__main__":
    main()
