#!/usr/bin/env python

'''
Uses the RNAStructure package 'bifold' on every combination of short sequences included in a <.fasta> file,
yeilding a <.csv> of the interaction energies between each combination at a given temperature.
'''

#Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import os
import subprocess
import glob
import itertools
import re

#Functions
def read_in_fasta(afasta):
    '''fasta file in, dictionary out'''
    fasta_sequences,fasta_dict = SeqIO.parse(open(afasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def populate_temp_directory(fasta_dict,directory):
    '''fill the directory with a fasta for every sequence'''
    os.mkdir(directory)
    os.chdir(directory)
    for name, sequence in fasta_dict.items():
        with open(name+'.fasta','w') as g:
            g.write('>'+name+'\n')
            g.write(sequence+'\n')
    os.chdir('..')

def generate_ct_files(directory,temp,DNA_flag):
    '''Walk back into the directory, get all perumations of fastas, dumb <.ct>s'''
    os.chdir(directory)
    all_fastas = itertools.combinations_with_replacement(sorted(glob.glob('*.fasta')),2)
    for fyle_1,fyle_2 in all_fastas:
        out_fyle = '~'.join([derp.split('.')[0] for derp in [fyle_1,fyle_2]])+'.ct'
        dna_bit = [] if DNA_flag == False else ['-d']
        command = ' '.join(['bifold',fyle_1,fyle_2,out_fyle,'-m','1','-t',temp]+dna_bit)#-m 1 because MFE.
        subprocess.call(command,shell=True)
    os.chdir('..')

def glean_ct_tops(directory):
    '''grab the top line of all files, return some sort of log'''
    energy_tank = {}
    os.chdir(directory)
    procced_info = subprocess.Popen('head -1 *.ct', shell=True, stdin=subprocess.PIPE, 
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True).stdout.read()
    read_data = [line for line in procced_info.split('\n') if line != '']
    for i in range(0,len(read_data),2):
        pair = read_data[i:i+2]
        key = pair[0][4:-4]
        try:
            value = re.findall("-\d+\.\d+", pair[1])[0]
        except IndexError:
            value = str(0)
        energy_tank[key] = value
    os.chdir('..')
    return energy_tank

def write_out_file(adict,outfyle):
    '''Takes information, writes it'''
    with open(outfyle,'w') as g:
        g.write(','.join(['File','Seq1','Seq2','Free_Energy'])+'\n')
        for k, v in adict.items():#I suppose we could produce a sorted file but w/e
            factor_1, factor_2 = k.split('~')[0],k.split('~')[1][:-3]
            if factor_1==factor_2:
                g.write(','.join([k,factor_1,factor_2,v])+'\n')
            else:
                g.write(','.join([k,factor_1,factor_2,v])+'\n')
                g.write(','.join([k,factor_2,factor_1,v])+'\n')

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def main():
    parser = argparse.ArgumentParser(description='Takes in a <.fasta> and uses Bifold on all sequence combinations')
    parser.add_argument('in_fasta',type=str,help='<.fasta> to pull sequences from')
    parser.add_argument('-name',type=str,default=None,help='Specify output name')
    parser.add_argument('-K',type=str, default='315.15', help='[default = 315.15] Kelvin temperature to fold at')
    parser.add_argument('-tempdir',type=str,default='temp',help='[default = temp] temporary directory to use')
    parser.add_argument('-savetemp',action='store_true',default=False,help = 'Do not delete the temporary directory when done')
    parser.add_argument('-DNA',action='store_true',default=False,help = 'Fold as DNA, not RNA')
    args = parser.parse_args()
    
    #Generate name
    nucl_tag = 'RNA' if args.DNA == False else 'DNA'
    s_name = '_'.join([args.in_fasta.split('.')[0],'bifolded',args.K,nucl_tag])+'.csv'
    outname = s_name if args.name == None else check_extension(args.name,'.csv')
    
    #Generate and write out.
    all_sequences = read_in_fasta(args.in_fasta)
    populate_temp_directory(all_sequences,args.tempdir)
    generate_ct_files(args.tempdir,args.K,args.DNA)
    information = glean_ct_tops(args.tempdir)
    write_out_file(information,outname)
    if args.savetemp == False:
        subprocess.call('rm -r temp',shell=True)


if __name__ == '__main__': 
    main()


