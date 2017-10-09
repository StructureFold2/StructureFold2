#!/usr/bin/python

'''
Returns max, average, standard deviation, and gini of reactivity for each transcript in <.react> files.
Operates on the entire directory at once, each set of columns will be named according to <.react> file names.
May have the output restricted by a coverage overlap list
Minumum length to operate on (flag m) is set to 10bp
Length to ignore from 3' end is set to 0 (flag n)
Gini function modified from https://planspacedotorg.wordpress.com/2013/06/21/how-to-calculate-gini-coefficient-from-raw-data-in-python/
The minimum length option works of the number of bases that are resolvable, i.e. at least N bases with a non NA value
'''

#Imports
import argparse
import glob
import numpy
from itertools import islice


#Functions

def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = None
    return good

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

def generate_reactivity_dictionary(react_fyles):
    '''Read in as many <.react> as you need to....'''
    all_data = {}
    #Read data in, builds a nested dictionary; missing keys are a problem for later
    for react_fyle in react_fyles:
        sub_data = read_derived_reactivities(react_fyle)
        for transcript, values in sub_data.items():
            all_data.setdefault(transcript, {})[react_fyle.replace('.react','')] = values
    return all_data

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def gini(list_of_values):
    '''Returns the Gini value of a list of values.'''
    try:
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area
    except ZeroDivisionError:
        return 'NA'

def generate_statistics(reactivity_info,react_fyles,trim_number,min_len):
    '''Generates the numbers for the different metrics by key, returns a new dictionary with the output'''
    order = sorted([x.replace('.react','') for x in react_fyles])
    ordered_statistics = {}
    for transcript, subdict in reactivity_info.items():
        new_value = []
        for item in order:
            try:
                numbers = [float(x) for x in subdict[item] if x!= 'NA'][:-trim_number] if trim_number else [float(x) for x in subdict[item] if x!= 'NA']
                if len(numbers) < min_len:
                    new_value.extend(['NA']*4)
                else:
                    try:
                        new_value.extend([max(numbers),numpy.average(numbers),numpy.std(numbers),gini(numbers)])
                    except ValueError:
                        new_value.extend(['NA']*4)
            except KeyError:
                new_value.extend(['NA']*4)
        ordered_statistics[transcript] = new_value
    return ordered_statistics

def write_out_file(outfyle,react_fyles,stats_dictionary):
    '''Write out the data in <.csv> format'''
    types,names = ['_max','_average','_std','_gini'],sorted([x.replace('.react','') for x in react_fyles])
    desc = [name+mod for name in names for mod in types]
    with open(outfyle,'w') as g:
        header = ['transcript']+desc
        g.write(','.join(header)+'\n')
        for k, v in stats_dictionary.items():
            fyle_lyne = ','.join([k]+[str(x) for x in v])
            g.write(fyle_lyne+'\n')

def main():
    parser = argparse.ArgumentParser(description='Generates a simple statistical report for all <.react> files in the directory.')
    parser.add_argument('-react',default = None, help='Operate on specific <.react> files', nargs='+')
    parser.add_argument('-restrict',default=None, help='Filter to these transcripts via coverage file')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-n',type=int, default=20, help='<int> [default = 20] ignore n last bp of reactivity',dest='trim')
    parser.add_argument('-m',type=int, default=10, help='<int> [default = 10] minimum length of transcript',dest='minlen')
    args = parser.parse_args()
    
    #Read in all queried <.react> files
    target_files = glob.glob('*.react') if args.react == None else args.react
    data = generate_reactivity_dictionary(target_files)
    
    #Read and apply restrictions if applicable
    restrict_dict = None if args.restrict == None else read_in_target_transcripts(args.restrict)
    if restrict_dict != None:
        filter_dictonary(data,restrict_dict)

    #Outfile nomenclature and assignment
    name_bits = sorted([x.replace('.react','') for x in target_files])
    ignored,minbases,suffix = [str(args.trim)+'trim'],[str(args.minlen)+'minlen'],['statistics.csv']
    default_name = '_'.join(name_bits+ignored+minbases+suffix)
    out_name = default_name if args.name == None else args.name
    
    #Generate Statistics, write out.
    out_data = generate_statistics(data,target_files,args.trim,args.minlen)
    write_out_file(out_name,target_files,out_data)

if __name__ == '__main__':
    main()
