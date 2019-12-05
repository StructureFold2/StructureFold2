#!/usr/bin/env python2

'''
Finds the 5' or 3' coverage of all transcripts in a .rtsc file.
It is highly recommended that only DMS(-) .rtsc files are used for input.

5' coverage is calculated by dividing the number of stops in the 5' most n number of nucleotides
(n is set by the user with -length) by the average number of stops for an equal sized region of the transcript.
To calculate 5' coverage, set FP as mode and set the value of length. -tp_l and -trim are not required.

3' coverage divides the number of stops in the 3' most n number of transcripts (n is again set by the user with -length)
by the number of stops in larger region of the 3' end of the transcript (set by the user with -tp_l).
It is recommended that 3' coverage is calculated following the trimming of increasing numbers of nucleotides from
the 3' end (set by the user with -trim) to assess the appropriate number of nucleotides to trim for all downstream analysis.
Any transcript that is shorter than the sum of the trim and tp_l values is returned as NA.
To calculate 3' coverage, set TP as mode and set the length, tp_l and trim values.

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
        if len(stops) + -parameters['trim'] + -parameters['tp_l'] > 0:
            try:
                trimmed_stops = stops[:-parameters['trim']]
                full_end = trimmed_stops[-parameters['tp_l']:]
                partial_end = trimmed_stops[-parameters['length']:]
                new_data[name] = average(partial_end)/average(full_end)
            except ZeroDivisionError:
                new_data[name] = 0.0
        else:
            new_data[name] = 'NA'
    extra_params = [str(parameters[q]) for q in ['length','tp_l','trim']]
    out_name = '_'.join([parameters['rtsc'].replace('.rtsc',''),'TP']+extra_params)
    out_name = check_extension(out_name,'.csv')
    return out_name, new_data

def write_coverage(coverage_dict,out_file,mode):
    '''Writes out in csv format'''
    header = ','.join(['transcript',mode+'_coverage'])
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
    write_coverage(coverage,out_file,args.mode)

if __name__ == '__main__':
    main()
