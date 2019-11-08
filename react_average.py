#!/usr/bin/env python

'''
This script averages any numbers of <.react> files.
It is imperative that all the files used to generate the <.react> used an identical reference in generation.
This module now impliclity checks for shared transcripts and lengths between the files before merging
and only outputs transcripts where everything checks out.
Module idea and original implementation by Joseph Waldron.
'''

#Imports
import argparse
from itertools import islice

#Functions
def read_in_react(react_fyle):
    '''Reads a <.react> file into a dictionary'''
    information = {}
    with open(react_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,stops = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x != 'NA' else 'NA' for x in stops.split()]
    print 'File: {} had a total of {} entries'.format(react_fyle,len(information))
    return information

def write_react(react_dictionary,sort_flag=False,outfile='data.react'):
    '''Takes a dictionary, writes to a file'''
    with open(outfile,'w') as g:
        if sort_flag == False:
            for transcript, value in react_dictionary.items():
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in value])+'\n')
        else:
            for transcript, value in sorted(react_dictionary.items()):
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in value])+'\n')

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def compress_lists(list_of_lists):
    '''Merges reactivity value lists'''
    out_list = []
    #We want to check if all the all the files have the same length associated with each key
    if len(set(map(len,list_of_lists))) == 1:
        for row in zip(*list_of_lists):
            #Note that we want to also check that all positions agree on numeric or NA
            if all(isinstance(item, str) for item in row):
                out_list.append('NA')
            #The data is easy to manipulate here as floats
            elif all(isinstance(item, float) for item in row):
                out_list.append(sum(row)/len(row))
            #But if there's some discrepancy where it isn't all NAs or all floats, oops!
            else:
                return None
        return out_list
    else:
        return None

def average_reacts(reacts):
    '''Reads in any number of reacts, runs compress_lists for the common keys'''
    new = {}
    all_reacts= map(read_in_react,reacts)
    #Generate common keys among all files, only work with these.
    common_keys = set.intersection(*map(set, all_reacts))
    for key in common_keys:
        new_values = compress_lists([sub[key] for sub in all_reacts])
        if new_values:
            new[key] = new_values
    return new

def main():
    parser = argparse.ArgumentParser(description='Combines <.react> files, typically replicates of the same sample')
    parser.add_argument('react',help='Input <.react> files', nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    args = parser.parse_args()
    
    #Generate name or assign the user provided name
    default_name = '_'.join(sorted([x.replace('.react','') for x in args.react]))+'_average.react'
    out_name = default_name if args.name == None else check_extension(args.name,'.react')
    
    #Average all the reacts
    averaged_reactivities = average_reacts(args.react)
    print 'Among all files, total entries merged: {}'.format(len(averaged_reactivities))
    
    #Write out the dictionary
    write_react(averaged_reactivities,args.sort,out_name)

if __name__ == '__main__': 
    main()
 