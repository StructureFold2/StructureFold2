#!/usr/bin/env python

'''
Take single files, or batch files, or both of coverage information from <.rtsc> files
and combine them into a single list of transcripts at or above a given threshold.
Note that same-named columns between files will be overridden, so do take care.
'''

#Imports
import glob
import os
import argparse

def read_coverage(coverage_file):
    ''''''
    all_the_things = {}
    with open(coverage_file,'r') as f:
        for line in f:
            if line.startswith('transcript'):
                file_keys = line.strip().split(',')[1:]
            else:
                line_components = line.strip().split(',')
                line_key,line_data = line_components[0],line_components[1:]
                assignments = zip(file_keys,line_data)
                for assignment in assignments:
                    all_the_things.setdefault(assignment[0],{})[line_key]=float(assignment[1])
    return all_the_things
            
def write_flat_list(master_dictionary,xthreshold=1):
    '''Writes a flat list, file is to be used by other scripts'''
    xoverlap,missingno = 0,0
    outname ='_'.join([group.replace('_coverage','') for group in sorted(master_dictionary.keys())]+['overlap',str(xthreshold)])+'.txt'
    common_keys = set([a for b in [sub.keys() for sub in master_dictionary.values()] for a in b])
    with open(outname,'w') as g:
        for key in common_keys:
            try:
                test_list = [(file_key,master_dictionary[file_key][key]) for file_key in sorted(master_dictionary.keys())]
                evaluated_test = [info[1]>=xthreshold for info in test_list]
                if all(evaluated_test):
                    g.write(key+'\n')
                    xoverlap+=1
            except KeyError:
                missingno +=1
    #Some sort of out metric
    print ''
    for filename,subdict in master_dictionary.items():
        print '{} had {} transcripts over {} coverage'.format(filename,len([x for x in subdict.values() if x >= xthreshold]),str(xthreshold))
    print 'Wrote file {} which had a total of {} mutually overlapped transcripts'.format(outname,str(xoverlap))
    print ''

def main():
    parser = argparse.ArgumentParser(description='Creates lists of transcripts with coverage overlaps above a certain threshold.')
    parser.add_argument('-n',type=float, default=1.0, help='[default = 1.0] coverage threshold to use',dest='threshold')
    parser.add_argument('-f', default = None, help='Specific <.csv>s to use', nargs='+',dest='csvs')
    parser.add_argument('-batchdir',action="store_true",default=False,help = 'Use all coverage <.csv> in the directory')
    args = parser.parse_args()
    #Initialize master coverage dictionary.
    all_coverages = {}
    #Get individual files
    if args.csvs != None:
        for fyle in args.csvs:
            sub_data = read_coverage(fyle)
            all_coverages.update(sub_data)
    #Get files in batch mode
    if args.batchdir == True:
        for fyle in sorted(glob.glob('*_coverage.csv')):
            sub_data = read_coverage(fyle)
            all_coverages.update(sub_data)
    
    #Write
    write_flat_list(all_coverages,args.threshold)

if __name__ == '__main__': 
    main()
 
