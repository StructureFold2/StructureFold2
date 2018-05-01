#!/usr/bin/env python

'''
Calculates comparative RMSD or NRMSD (normalized) between two <.react> files.
'''

import argparse
from itertools import islice

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def read_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!= 'NA' else 'NA' for x in reactivities.split()]
    return information

def calculate_RMSD(vector_1,vector_2):
    '''Calculates Standard RMSD on two vectors of numbers of the same length'''
    #Check to see the vectors are of equal length.
    if len(vector_1)!=len(vector_2):
        return 'NA'
    eq_top,total_n = 0,0
    #Sum the squares of the differences/total N where possible.
    for matched in zip(vector_1,vector_2):
        try:
           difference = (float(matched[0])-float(matched[1]))**2
           eq_top+=difference
           total_n+=1
        except:
            continue
    #Divide total difference by total N, take square root.
    result = (eq_top/total_n)**0.5 if total_n != 0 else 'NA'
    return result

def calculate_NRMSD(vector_1,vector_2):
    '''Calculates Standard RMSD on two vectors of numbers of the same length'''
    #Check to see the vectors are of equal length.
    if len(vector_1)!=len(vector_2):
        return 'NA'
    eq_top,total_n,values = 0,0,[]
    #Sum the squares of the differences and total N where possible, sum all numbs.
    for matched in zip(vector_1,vector_2):
        try:
           difference = (float(matched[0])-float(matched[1]))**2
           eq_top+=difference
           total_n+=1
           values.extend([float(matched[0]),float(matched[1])])
        except:
            continue
    #Divide total difference by total N, take square root, apply normalization.
    average_react = sum(values)/2/total_n if total_n !=0 else 0
    result = (eq_top/total_n)**0.5/average_react if total_n != 0 and average_react !=0 else 'NA'
    return result

def mass_calculate_RMSD(data1,data2):
    '''Runs the RMSD funciton on two dictionaries of mutual vectors'''
    values = {}
    iteration_keys = sorted(list(set(data1.keys()).intersection(set(data2.keys()))))
    for key in iteration_keys:
        values[key] = calculate_RMSD(data1[key],data2[key])
    return values

def mass_calculate_NRMSD(data1,data2):
    '''Runs the NRMSD funciton on two dictionaries of mutual vectors'''
    values = {}
    iteration_keys = sorted(list(set(data1.keys()).intersection(set(data2.keys()))))
    for key in iteration_keys:
        values[key] = calculate_NRMSD(data1[key],data2[key])
    return values

def write_out_results(info_dict,outfyle,column_name):
    '''Writes out the values as a <.csv>'''
    with open(outfyle,'w') as g:
        g.write(','.join(['transcript',column_name])+'\n')
        for transcript, value in info_dict.items():
            g.write(','.join([transcript,str(value)])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Calculates RMSD/NRMSD between two <.react> files')
    parser.add_argument('react1',help='<.react> file 1')
    parser.add_argument('react2',help='<.react> file 2')
    parser.add_argument('-normalize',action='store_true',default=False,help = 'Normalize output RMSD values')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-restrict',default=None, help='Restrict output to these transcripts')
    args = parser.parse_args()
    
    #Create output name
    end_suffix = ['rmsd.csv'] if args.normalize == False else ['nrmsd.csv']
    base_name = [x.replace('.react','') for x in [args.react1,args.react2]]
    out_name = '_'.join(base_name+end_suffix) if args.name == None else args.name
    
    #Read data in
    control,experimental = read_reactivities(args.react1),read_reactivities(args.react2)
    
    #Truncate seqs to those with good coverage in both conditions if user provides a list
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(control,covered)
        filter_dictonary(experimental,covered)
    
    #Process data
    data = mass_calculate_RMSD(control,experimental) if args.normalize == False else mass_calculate_NRMSD(control,experimental)
    
    #Write out data
    col_name = 'rmsd' if args.normalize == False else 'nrmsd'
    write_out_results(data,out_name,col_name)

if __name__ == '__main__':
    main()
