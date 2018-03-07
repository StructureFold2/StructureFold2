#!/usr/bin/env python

'''
Applies random.shuffle() to the sequences of a <.fasta> file, generating a
somewhat random prior for each sequence with the same length and nucleotide composition. 
Thus, these sequences can be folded, and the predicted structure and MFE of the random sequence
compared to the native sequence, where any difference may point to the sequence itself being 
potentially important rather than just length and composition
'''

#Imports
import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Functions
def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def write_out_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for name,seq in sorted(info.items()):
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n')

def shuffler(fasta_dict,pfix):
    '''shuffles'''
    new_dict = {}
    for name,sequence in fasta_dict.items():
        xlist = list(sequence)
        random.shuffle(xlist)
        name = name if pfix == None else '_'.join([pfix,name])
        new_dict[name] = ''.join(xlist)
    return new_dict

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Shuffles <fasta> sequences')
    parser.add_argument('fasta',type=str,help='<.fasta> file')
    parser.add_argument('-name',default='outfile.fasta', help='Specify output file name')
    parser.add_argument('-prefix',type=str,default=None, help='Prefix randomized seq names')
    args = parser.parse_args()
    
    #Read in fasta
    sequences = read_fasta(args.fasta)
    
    #Create new fasta dictionary
    new_sequences = shuffler(sequences,args.prefix)
    
    #Write new dictionary
    write_out_fasta(new_sequences,args.name)

if __name__ == '__main__': 
    main()
