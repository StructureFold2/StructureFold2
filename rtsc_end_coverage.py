#!/usr/bin/env python2

'''
Finds the 5' or 3' coverage of items in a .rtsc file
Idea and some code by Joseph Waldron
'''

#Imports
import argparse
import itertools

#Functions
def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary, transcript_name:[list of stop numbers]'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(itertools.islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    print '{} yielded {} entries'.format(rtsc_fyle,len(information))
    return information

def average(alist):
    '''Empty Docstring'''
    return sum(alist)/float(len(alist))

def check_extension(astr,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astr if astr.endswith(extension) else astr+extension
    return out_string

def five_prime_coverage(rtsc_data,parameters):
    '''Calculates FP coverage by Joseph's reckoning'''
    new_data = {}
    for name, stops in rtsc_data.items():
        try:
            new_data[name] = sum(stops[:parameters['length']])/(average(stops)*parameters['length'])
        except ZeroDivisionError:
            new_data[name] = 0.0
    extra_params = [str(parameters[q]) for q in ['length']]
    out_name = '_'.join([parameters['rtsc'].replace('.rtsc',''),'FP']+extra_params)
    out_name = check_extension(out_name,'.csv')
    return out_name,new_data

def three_prime_coverage(rtsc_data,parameters):
    '''Calculates TP coverage by Joseph's reckoning'''
    new_data = {}
    for name, stops in rtsc_data.items():
        try:
            trimmed_stops = stops[:-parameters['trim']]
            full_end = trimmed_stops[-parameters['tp_l']:]
            partial_end = trimmed_stops[-parameters['length']:]
            new_data[name] = average(partial_end)/average(full_end)
        except ZeroDivisionError:
            new_data[name] = 0.0
    extra_params = [str(parameters[q]) for q in ['length','tp_l','trim']]
    out_name = '_'.join([parameters['rtsc'].replace('.rtsc',''),'TP']+extra_params)
    out_name = check_extension(out_name,'.csv')
    return out_name, new_data

def write_coverage(coverage_dict,out_file):
    '''Writes out in csv format'''
    header = ','.join(['transcript','coverage'])
    with open(out_file,'w') as g:
        g.write(header+'\n')
        for transcript, value in sorted(coverage_dict.items(), key=lambda x: x[1],reverse=True):
            g.write(','.join([transcript,str(value)])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates a <.csv> of 5\' or 3\' end coverage from an <.rtsc> file')
    parser.add_argument('rtsc',help='file to operate on')
    parser.add_argument('mode',type=str.upper,help='5\' or 3\' Prime',choices = ['FP','TP'])
    parser.add_argument('-length',type=int,default=50,help='[5\'& 3\', default=50] Number of bases from the end')
    parser.add_argument('-tp_l',type=int,default=300,help='[3\', default=300] Number of comparative bases')
    parser.add_argument('-trim',type=int,default=30,help='[3\', default=30] Number of bases to ignore')
    parser.add_argument('-name',type=str,default=None,help='outfile name, overrides auto-gen name')
    args = parser.parse_args()
    functs = {'FP':five_prime_coverage,'TP':three_prime_coverage}
    
    #Read in an rtsc.
    rtsc = read_in_rtsc(args.rtsc)
    
    #Calculate coverage on the given end.
    out_file, coverage = functs[args.mode](rtsc,vars(args))
    
    #Choose final out name
    out_file = out_file if args.name == None else check_extension(args.name,'.csv')
    
    #Write Out
    write_coverage(coverage,out_file)

if __name__ == '__main__':
    main()