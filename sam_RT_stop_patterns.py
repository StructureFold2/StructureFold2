#!/usr/bin/env python
#David Tack

'''
Collects and reports the frequency of all di-nucleotide motifs in a <.sam> file. If the corresponding <.fasta>
file is included, the motif of the actual read stop is collected separately and also reported, such that it can be
compared to the overall transcriptome 'prior' for any bias.
'''

#Imports
import argparse
import sys
import glob
from collections import Counter
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Functions
def sequence_to_motif_dist(sequence,window=2,step=1):
    '''Returns a counter of the sequence motifs of a window of n and a step of n'''
    value_windows = Counter([sequence[i:i+window] for i in xrange(0, len(sequence)-(window-1), step)])
    print value_windows

def sequence_dinucleotide_freq(seq):
    '''Simple dinucelotide'''
    return Counter([seq[n-1:n+1] for n in range(1,len(seq))])

def read_fasta(genome_fasta):
    '''Reads a <.fasta> file, returns a dictionary.'''
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def clean_round(afloat):
    '''takes a float, reuturn a clean 2 decimal point string'''
    a,b= str(afloat).split('.')
    return '.'.join([a,b[0:2]])

def generate_prior(sam_file):
    '''Takes a <.sam> file, returns a counter of distributions'''
    motifs = Counter()
    with open(sam_file,'r') as f:
        for line in f:
            bits = line.strip().split()
            my_counter = sequence_dinucleotide_freq(bits[9])
            motifs = motifs+my_counter
    return motifs

def batch_prior(alist):
    '''batches up generate_prior'''
    batch_priors = {}
    for sam_fyle in alist:
        batch_priors[sam_fyle] = generate_prior(sam_fyle)
    return batch_priors

def batch_all(alist,sequence_reference):
    '''batches both'''
    batch_priors,batch_hits = {},{}
    for sam_file in alist:
        prior,hits = generate_distributions(sam_file,sequence_reference)
        batch_priors[sam_file] = prior
        batch_hits[sam_file] = hits
    return batch_priors, batch_hits

def generate_distributions(sam_file,sequence_dict):
    '''Generate both disributions'''
    motifs,splits = Counter(),Counter()
    with open(sam_file,'r') as f:
        for line in f:
            bits = line.strip().split('\t')
            gene_sequence,position,my_counter = sequence_dict[bits[2]],int(bits[3]),sequence_dinucleotide_freq(bits[9])
            motifs = motifs+my_counter
            if position != 1:
                splits[gene_sequence[position-2:position]]+=1
    return motifs,splits

def write_out_dist(prior_counter,out_file,report_threshold,source_file,d_counter=None):
    '''Writes out a twin distribution potentially with di-nucleotide RT Stop info'''
    base_header,extra_header = ['file','motif','prior_percent','prior_count'],['stop_percent','stop_count']
    header = base_header + extra_header if d_counter else base_header
    with open(out_file,'w') as g:
        g.write(','.join(header)+'\n')
        for k, v in sorted(prior_counter.items()):
            if v > report_threshold:
                part_1 = [source_file,k,clean_round(float(v)/sum(prior_counter.values())*100),str(v)]
                part_2 = [clean_round(float(d_counter[k])/sum(d_counter.values())*100),str(d_counter[k])] if d_counter else []
                g.write(','.join(part_1+part_2)+'\n')

def write_out_multi(prior_dict,out_file,report_threshold,stop_dict=None):
    '''Writes any number of distributions of priors or priors and stops'''
    base_header,extra_header = ['file','motif','prior_percent','prior_count'],['stop_percent','stop_count']
    header = base_header + extra_header if stop_dict else base_header
    with open(out_file,'w') as g:
        g.write(','.join(header)+'\n')
        for source,subdict in sorted(prior_dict.items()):
            for k, v in sorted(subdict.items()):
                if v > report_threshold:
                    part_1 = [source,k,clean_round(float(v)/sum(subdict.values())*100),str(v)]
                    part_2 = [clean_round(float(stop_dict[source][k])/sum(stop_dict[source].values())*100),str(stop_dict[source][k])] if stop_dict else []
                    g.write(','.join(part_1+part_2)+'\n')

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string


#Main Function
def main():
    parser = argparse.ArgumentParser(description='Analyses Di-nucleotide motif patterns from <.sam> files')
    parser.add_argument('-single',default = None, help = 'Operate on this single file, rather than the directory')
    parser.add_argument('-index',type=str,default = None, help='<.fasta> file to generate RT stop motifs')
    parser.add_argument('-min',type=int,default = 1, help='minimum occurances of a motif to report',dest='min_m')
    parser.add_argument('-name',type=str,default = 'outfile.csv', help='[default = outfile.csv] Name of the output')
    args = parser.parse_args()
    
    out_name = check_extension(args.name,'.csv')
    
    #Batch mode
    if args.single == None:
        sam_files = glob.glob('*.sam')

        if args.index:
            sequences = read_fasta(args.index)
            batch_priors,batch_stops = batch_all(sam_files,sequences)
            write_out_multi(batch_priors,out_name,args.min_m,batch_stops)

        else:
            batch_priors = batch_prior(sam_files)
            write_out_multi(batch_priors,out_name,args.min_m)
            
    #Single Mode
    else:
        if args.index:
            sequences = read_fasta(args.index)
            single_motifs,single_stops = generate_distributions(args.single, sequences)
            write_out_dist(single_motifs,out_name,args.min_m,args.single,single_stops)
        
        else:
            single_motifs = generate_prior(args.single)
            write_out_dist(single_motifs,out_name,args.min_m,args.single)


if __name__ == '__main__': 
    main()
