#!/usr/bin/python

'''
Extracts a single transcript from a <.react> and a <.fasta> file, reformats it into an easy to read <.csv>
'''

#Imports
import argparse
import os
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Functions
def read_derived_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            reactivities = [float(x) if x != 'NA' else 'NA' for x in reactivities.split()]
            information[transcript] = reactivities
    return information

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def write_out_csv(sequence,values,outfile='stats.csv'):
    '''Writes out the data'''
    with open(outfile,'w') as g:
        g.write(','.join(['Position','Nucleotide','Reactivity'])+'\n')
        lines = zip(range(1,len(sequence)+1),sequence,values)
        for line in lines:
            g.write(','.join([str(x) for x in line])+'\n')


def batch_write_out_csv(sequence_dict,values_dict,name_suffix):
    '''batch writes a lot of csv'''
    q_keys = set(sequence_dict.keys()).intersection(set(values_dict.keys()))
    for key in q_keys:
        new_fyle = '_'.join([key,name_suffix])+'.csv'
        sequence, values = sequence_dict[key],values_dict[key]
        if len(sequence) == len(values):
            with open(new_fyle,'w') as g:
                g.write(','.join(['Position','Nucleotide','Reactivity'])+'\n')
                lines = zip(range(1,len(sequence)+1),sequence,values)
                for line in lines:
                    g.write(','.join([str(x) for x in line])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Reformats reactivity values to <.csv> format')
    parser.add_argument('react',type=str,help='<.react> file to pull values from')
    parser.add_argument('fasta',type=str,help='<.fasta> file used to generate the <.react>')
    parser.add_argument('-transcript',default=None,type=str,help='Specific transcript to reformat')
    parser.add_argument('-all_transcripts',action='store_true',help='Reformat all, output to new directory')
    args = parser.parse_args()
    
    #Read in files to query
    sequences, reactivities = read_in_fasta(args.fasta),read_derived_reactivities(args.react)
    
    if args.transcript and not args.all_transcripts:
        #Check Values
        seq = sequences[args.transcript] if args.transcript in sequences else None
        reacts = reactivities[args.transcript] if args.transcript in reactivities else None
        #Auto-generate name
        outname= '_'.join([args.transcript,args.react.replace('.react','')])+'.csv'
        #Write
        if seq != None and reacts !=None:
            if len(seq) == len(reacts):
                write_out_csv(seq,reacts,outname)
            else:
                print 'Sequence length does not match number of reactivities for {}'.format(args.transcript)
        else:
            print 'One or more of the files did not contain entry {}'.format(args.transcript)
    
    elif args.all_transcripts and not args.transcript:
        #Directory
        new_dir = args.react.replace('.react','')+'_'+'all_csvs'
        check= os.listdir('.')
        if new_dir not in check:
            os.mkdir(new_dir)
            os.chdir(new_dir)
            batch_write_out_csv(sequences,reactivities,args.react.replace('.react',''))
            os.chdir('..')
        else:
            print '{} folder already exists. Remove or rename and try again.'.format(new_dir)
    
    else:
        print 'Check options, use only one of -all_transcripts and -transcript'

if __name__ == '__main__':
    main()
