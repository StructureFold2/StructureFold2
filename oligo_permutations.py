#!/usr/bin/env python

#David Tack
'''Takes a base oligo or generates one. Chains mutations in sequentially, have any number of forks'''

#Imports
import argparse
import random

#Functions
def generate_random_seq(alpha,l=20):
    '''Generates a random nucleotide sequence'''
    return ''.join([random.choice(alpha) for z in range(0,l)])

def add_snp(seq,choices):
    '''Introduces a single SNP, pulls from choices'''
    temp,index = list(seq),random.randint(0,len(seq)-1)
    before =temp[index]
    new = random.choice([item for item in choices if item != before])
    temp[index] = new
    return ''.join(temp)

def chain_iterate(sequence,stepz,alphabet,I=1):
    '''does an n-step chain of iterations, no backsteps'''
    new_seqs = {'delirium':sequence}
    frozen_save = sequence
    while len(new_seqs) < stepz+1:
        cold_save = add_snp(frozen_save,alphabet)
        if cold_save not in new_seqs.values():
            new_seqs['_'.join(['chain',str(I).zfill(4)])] = cold_save
            frozen_save = cold_save
            I+=1
    del new_seqs['delirium']
    return new_seqs

def add_snp_restrict(seq,banned,choices):
    '''Introduces a single SNP, pulls from choices'''
    temp = list(seq)
    index = random.randint(0,len(seq)-1)
    while index in banned:
        index = random.randint(0,len(seq)-1)
    before = temp[index]
    new = random.choice([item for item in choices if item != before])
    temp[index] = new
    return ''.join(temp)

def chain_iterate_restrict(sequence,stepz,nonobases,alphabet,I=1):
    '''does an n-step chain of iterations, no backsteps'''
    new_seqs = {'delirium':sequence}
    frozen_save = sequence
    while len(new_seqs) < stepz+1:
        cold_save = add_snp_restrict(frozen_save,nonobases,alphabet)
        if cold_save not in new_seqs.values():
            new_seqs['_'.join(['chain',str(I).zfill(4)])] = cold_save
            frozen_save = cold_save
            I+=1
    del new_seqs['delirium']
    return new_seqs

def run_forks(sequence,n_forks,n_chain,alpha,restrict_bases=None):
    '''Runs the Forks'''
    master = {}
    for i in range(0,n_forks):
        sub = chain_iterate(sequence,n_chain,alpha) if restrict_bases == None else chain_iterate_restrict(sequence,n_chain,restrict_bases,alpha) 
        for key, value in sorted(sub.items()):
            if value not in master.items():
                zkey = '_'.join(['fork',str(i+1).zfill(4),key])
                master[zkey] = value
            else:
                continue
    master['original'] = sequence
    return master

def write_out_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for name,seq in sorted(info.items()):
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Creates <.fasta> of sequence variants.')
    parser.add_argument('-prior',type=str,default=None,help='[default = random] Nucleotide Sequence to iterate on')
    parser.add_argument('-length',type=int,default=20, help='[default = 20] Number of bases for random prior')
    parser.add_argument('-chain',type=int,default=20, help='[default = 20] Number of iterations from base seq')
    parser.add_argument('-fork',type=int,default=8, help='[default = 8] Number of forks from base seq')
    parser.add_argument('-name',type=str,default='out.fa', help='[default = out.fa] Name the output')
    parser.add_argument('-alphabet',type=str,default='ACGT',help='[default = ACGT] Alphabet to use')
    parser.add_argument('-r',type=int,default = None, help='Bases numbers which may not be permuted', nargs='+',dest='same')
    args = parser.parse_args()
    
    #Generate sequence if none provided
    seed = generate_random_seq(args.alphabet,args.length) if args.prior == None else args.prior.upper()
    
    #Do it
    fasta_dict = run_forks(seed,args.fork,args.chain,args.alphabet,args.same)
    
    #Write OUT
    write_out_fasta(fasta_dict,args.name)


if __name__ == '__main__': 
    main()
