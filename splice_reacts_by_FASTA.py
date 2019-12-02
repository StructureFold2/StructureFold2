#!/usr/bin/env python2

'''
This script will splice a react file into three react files, one each for the 5'UTR, the CDS and the 3'UTR.
It requires an already spliced FASTA and will splice the react file exactly as the spliced FASTA files.
For transcripts that have a CDS but not a 5'UTR or 3'UTR, it will still add the spliced CDS reactivities to the CDS file,
but these transcripts will not be in the 5'UTR or 3'UTR file.
'''

from itertools import islice
from Bio import SeqIO
import argparse

def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(afasta),'fasta')
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq).upper()
    return fasta_dict

def annotate_UTRs(fp_FASTA, CDS_FASTA, tp_FASTA):
    '''takes the FASTA for each spliced region and creates a dictionary
    with transcript ID as key and a tuple containing 5'UTR, CDS and 3'UTR lengths as the value'''
    adict = {}
    for k,v in CDS_FASTA.items():
        CDS_length = len(v)
        try:
            fpUTR_length = len(fp_FASTA[k])
        except KeyError:
            fpUTR_length = 0
        try:
            tpUTR_length = len(tp_FASTA[k])
        except KeyError:
            tpUTR_length = 0
        adict[k] = (fpUTR_length, CDS_length, tpUTR_length)
    return adict

def write_spliced_react_files(reactivity_file, annotation_dict):
    '''takes a react file, splices it into 5'UTRs/CDS/3'UTR and writes a new react file for each region'''
    FP_OUT,CDS_OUT,TP_OUT = [reactivity_file.split('.')[0] + n for n in ['_fpUTR.react','_CDS.react','_tpUTR.react']]
    with open(reactivity_file,'r') as f,open(FP_OUT,'w') as a,open(CDS_OUT,'w') as b,open(TP_OUT,'w') as c:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            expanded_reactivities = reactivities.split('\t')
            if transcript in annotation_dict:
                fp_length = annotation_dict[transcript][0]
                CDS_length = annotation_dict[transcript][1]
                tp_length = annotation_dict[transcript][2]
                if len(expanded_reactivities) == fp_length + CDS_length + tp_length:
                    if fp_length > 0:
                        FP_UTR = expanded_reactivities[:fp_length]
                        a.write(transcript+'\n')
                        a.write('\t'.join(FP_UTR)+'\n')
                    CDS = expanded_reactivities[fp_length:-tp_length]
                    b.write(transcript+'\n')
                    b.write('\t'.join(CDS)+'\n')
                    if tp_length > 0:
                        TP_UTR = expanded_reactivities[-tp_length:]
                        c.write(transcript+'\n')
                        c.write('\t'.join(TP_UTR)+'\n')
                else:
                    print transcript + " was not spliced as the lengths of the spliced FASTA files did not equal the length of unspliced reactivities"

def main():
    parser = argparse.ArgumentParser(description = 'Splices REACT files into regions based on on the lengths of the sequences in the spliced FASTA files')
    parser.add_argument('react', type = str, help = 'REACT file')
    parser.add_argument('fp_FASTA', type = str, help = '5\'UTR FASTA file')
    parser.add_argument('CDS_FASTA', type = str, help = 'CDS FASTA file')
    parser.add_argument('tp_FASTA', type = str, help = '3\'UTR FASTA file')
    args = parser.parse_args()
    
    fpUTRs, CDSs, tpUTRs = read_in_fasta(args.fp_FASTA), read_in_fasta(args.CDS_FASTA), read_in_fasta(args.tp_FASTA)
    annotated_transcripts = annotate_UTRs(fpUTRs, CDSs, tpUTRs)
    write_spliced_react_files(args.react, annotated_transcripts)


if __name__ == '__main__': 
    main()