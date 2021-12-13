#!/usr/bin/env python2

'''
By default, this will compare every <.rtsc> in the directory to the -index fasta file, and generate a small report.
A custom set of files can instead be given to the program.
'''

#Imports
import glob
import argparse
from collections import Counter
from sf2libs.structure_io import read_fasta, read_rtsc, check_extension

#Functions
def rtsc_to_specificity(rtsc_fyle,fasta_index):
    '''Returns a Counter object of the frequencies of the nucleotide specificities'''
    Xcounter = Counter()
    rtsc_dictionary = read_rtsc(rtsc_fyle)
    for transcript, stops in rtsc_dictionary.items():
        seq = fasta_index[transcript].upper()[:-1]
        matched_values = zip(list(seq),stops[1:])
        for pair in matched_values:
            Xcounter[pair[0]]+=int(pair[1])
    return Xcounter
            
def write_specificity_data(info,outfyle='specificity.csv'):
    '''Dumps to a .csv'''
    header_keys = sorted([x.replace('.rtsc','') for x in info.keys()])
    types = ['_count','_specificity']
    header = ','.join(['base']+[fyle+mod for fyle in header_keys for mod in types])
    f_keys = sorted(info.keys())
    A = ['A']+[(info[q]['A'],info[q]['A']/float(sum(info[q].values()))) for q in f_keys]
    T = ['T']+[(info[q]['T'],info[q]['T']/float(sum(info[q].values()))) for q in f_keys]
    G = ['G']+[(info[q]['G'],info[q]['G']/float(sum(info[q].values()))) for q in f_keys]
    C = ['C']+[(info[q]['C'],info[q]['C']/float(sum(info[q].values()))) for q in f_keys]
    x_data = [A,T,G,C]
    with open(outfyle,'w') as g:
        g.write(header+'\n')
        for data in x_data:
            g.write(','.join(str(y) for x in data for y in x)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Analyzes native/reagent nucleotide stop specificity')
    parser.add_argument('-index',type=str,help='<.fasta> file used to generate the <.rtsc>')
    parser.add_argument('-rtsc',default = None, help='Operate on specific <.rtsc>', nargs='+')
    parser.add_argument('-name',default = None, help='Specify output file name')
    args = parser.parse_args()
    
    #Read in fasta sequences
    seqs = read_fasta(args.index)

    #Batch Mode
    if args.rtsc == None:
        #Generate Out Name
        fyles,info = sorted(glob.glob('*.rtsc')),{}
        default_name = '_'.join([x.replace('.rtsc','') for x in fyles])+'_specificity.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.rtsc')
        for fyle in fyles:
            info[fyle] = rtsc_to_specificity(fyle,seqs)
        write_specificity_data(info,out_name)
    
    #Specific Mode
    if args.rtsc != None:
        #Generate Out Name
        info = {}
        default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'_specificity.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.rtsc')
        for fyle in args.rtsc:
            info[fyle] = rtsc_to_specificity(fyle,input_fasta)
        write_specificity_data(info,out_name)

if __name__ == '__main__': 
    main()
