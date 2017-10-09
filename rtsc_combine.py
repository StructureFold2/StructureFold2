#!/usr/bin/env python

'''
This script combines any numbers of <.rtsc> files.
It is imperative that all the files used to generate the <.rtsc> used an identical reference in generation.

'''

#Imports
import argparse
from itertools import islice

#Functions
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

def write_rtsc(rtsc_dictionary,sort_flag=False,outfile='data.rtsc'):
    '''Takes a dictionary, writes to a file'''
    with open(outfile,'w') as g:
        if sort_flag == False:
            for transcript, value in rtsc_dictionary.items():
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in value])+'\n\n')
        else:
            for transcript, value in sorted(rtsc_dictionary.items()):
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in value])+'\n\n')

def merge_rtsc(rtsc_lyst):
    '''Takes a list of <.rtsc> files, returns a dictionary that is sum of the RT stops'''
    master_dictionary= {}
    for rtsc in sorted(rtsc_lyst):
        data = read_in_rtsc(rtsc)
        for transcript, values in data.items():
            if transcript in master_dictionary:
                current_values = master_dictionary[transcript]
                master_dictionary[transcript] = [sum(spot) for spot in zip(current_values,values)]
            else:
                master_dictionary[transcript] = values
    return master_dictionary

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def main():
    parser = argparse.ArgumentParser(description='Combines <.rtsc> files, typically replicates of the same sample')
    parser.add_argument('rtsc',help='Input <.rtsc> files', nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    args = parser.parse_args()
    #Generate name or assign the user provided name
    default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'.rtsc'
    out_name = default_name if args.name == None else check_extension(args.name,'.rtsc')
    #Pool all <.rtsc> into a dictionary
    all_stops = merge_rtsc(args.rtsc)
    #Write out the dictionary
    write_rtsc(all_stops,args.sort,out_name)


if __name__ == '__main__': 
    main()
 
 
