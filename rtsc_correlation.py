#!/usr/bin/env python

'''
Reformats <.rtsc> files such that correlation may be easily caluclated, either genome-wide or on a per transcript basis.
Output is a <.csv>.
You may filter output by using a coverage overlap. where you only want to get the correlation between transcripts mutually above a certain threshold.
You may filter with a nucleotide specificity filter, but this requires use of the <.fasta> mapped against via Bowtie.
'''

#Imports
import argparse
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Functions
def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = 'NULL'
    return good

def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information

def generate_repeatabilty_dictionary(rtsc_fyles):
    '''Read in as many <.rtsc> as you need to....'''
    all_data = {}
    #Read data in, builds a nested dictionary
    for rtsc_fyle in rtsc_fyles:
        sub_data = read_in_rtsc(rtsc_fyle)
        for transcript, values in sub_data.items():
            all_data.setdefault(transcript, {})[rtsc_fyle.replace('.rtsc','')] = values
    return all_data
            
def write_out_repeatability(stops,in_files,outfile='OUT.csv'):
    '''Writes out a nested RTSC dictionary'''
    with open(outfile,'w') as g:
        header = ','.join(['transcript','position']+sorted([x.replace('.rtsc','') for x in in_files]))+'\n'
        g.write(header)
        for transcript, subdict in stops.items():
            all_data,pos = zip(*[subdict[key][1:] for key in sorted(subdict.keys())]),1
            for item in all_data:
                g.write(','.join([transcript,str(pos)]+[str(x) for x in item])+'\n')
                pos+=1

def write_out_repeatability_spec(stops,in_files,seqs,specificity,outfile='OUT.csv'):
    '''Writes out a nested RTSC dictionary with a given specificity'''
    accepted = list(specificity)
    with open(outfile,'w') as g:
        header = ','.join(['transcript','position','base']+sorted([x.replace('.rtsc','') for x in in_files]))+'\n'
        g.write(header)
        for transcript, subdict in stops.items():
            base_seq = list(seqs[transcript])[:-1]
            x_frame = [subdict[key][1:] for key in sorted(subdict.keys())]
            g_frame = [base_seq]+x_frame
            all_data,pos = zip(*g_frame),1
            for item in all_data:
                if item[0] in accepted:
                    g.write(','.join([transcript,str(pos)]+[str(x) for x in item])+'\n')
                pos+=1

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(afasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def main():
    parser = argparse.ArgumentParser(description='Reformats <.rtsc> for easy correlation analysis')
    parser.add_argument('rtsc',help='Input <.rtsc> files', nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-fasta',default=None, help='<.fasta> to apply specificity')
    parser.add_argument('-spec',default=None, help='[ACGT] Nucleotide Specifictiy')
    parser.add_argument('-restrict',default=None, help='Filter to these transcripts via coverage file')
    args = parser.parse_args()
    
    #Outfile nomenclature
    default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'_correlation.csv'
    if args.spec != None:
        new_suffix = '_'+''.join(sorted(list(args.spec)))+'spec_correlation.csv'
        default_name = default_name.replace('_correlation.csv',new_suffix)
    out_name = default_name if args.name == None else args.name
    
    #Restrictions, if they exist
    restrict_dict = None if args.restrict == None else read_in_target_transcripts(args.restrict)
    
    #Read in data, filter out transcripts if applicable
    data = generate_repeatabilty_dictionary(args.rtsc)
    if restrict_dict != None:
        filter_dictonary(data,restrict_dict)

    #No specificity path
    if args.fasta == None and args.spec == None:
        write_out_repeatability(data,args.rtsc,out_name)
    
    #Specificity path
    elif args.fasta != None and args.spec != None:
        sequences = read_in_fasta(args.fasta)
        write_out_repeatability_spec(data,args.rtsc,sequences,args.spec,out_name)
    
    else:
        print 'Invalid command combination.'
        print '-spec and -fasta must be invoked together'
        print ''


if __name__ == '__main__': 
    main()
 
