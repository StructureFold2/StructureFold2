#!/usr/bin/python
#David Tack

'''
Performs various binning functions on <.react> files. 
The number and minimum enforced bin size are modifiable. Three Modes:
Single Mode pulls one targeted transcript from the <.react>
Multi Mode takes a flat <.txt> list of transcripts, and bins all it can together, taking the average of averages at each bin.
All mode takes all transcripts in a <.react> file that pass bin criteria.
'''

#Imports
import argparse
import numpy
from itertools import islice

#Functions
def read_in_derived_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = reactivities.split('\t')
    return information

def split_values(alyst, n):
    '''Returns n mostly even chunks of a list'''
    k, m = len(alyst) / n, len(alyst) % n
    return [alyst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)]

def average_sub_list(alyst):
    '''returns the average of each sublist'''
    for item in alyst:
        sub = [float(x) for x in alyst if x != 'NA']
        if sub == []:
            return 0.0
        else:
            return sum(sub)/len(sub)

def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = None
    return good

def dump_list_as_file(alist,source_file,transcript):
    '''Write out bins'''
    outfyle = source_file.split('.')[0]+'_'+transcript+'_binned.csv'
    with open(outfyle,'w') as f:
        f.write('value\n')
        for item in alist:
            f.write(str(item)+'\n')

def add_bins(lyst_a,lyst_b):
    '''Sums two flat lists'''
    return [sum(x) for x in zip(lyst_a,lyst_b)]
    
def assert_bin_sizes(lyst_of_lysts,threshold_size):
    '''Returns true or false depending if every sub list contains at least n values'''
    return all([len(x) >= threshold_size for x in lyst_of_lysts])

def main():
    parser = argparse.ArgumentParser(description='Generates reactivity bins of transcripts from a <.react> file')
    parser.add_argument('react', type=str, help='<.react> file to pull values from')
    parser.add_argument('-bins',type=int, default=100, help='[default = 100] bins to create')
    parser.add_argument('-size',type=int, default=10, help='[default = 10] minimum number of values per bin')
    parser.add_argument('-single', type=str,default=None, help='Target transcript to bin from the file')
    parser.add_argument('-multi', type=str,default=None, help='Flat <.txt> of transcripts to bin, one per line')
    parser.add_argument('-all_transcripts',action='store_true', help='Amalgamate all transcripts in the <.react> file')
    args = parser.parse_args()
    
    #Read in the <.react> file
    transcripts = read_in_derived_reactivities(args.react)
    
    #Single Mode
    if args.single and not args.all_transcripts and not args.multi:
        try:
            vector = split_values(transcripts[args.single],args.bins)
            if assert_bin_sizes(vector,args.size):
                vector_average = [average_sub_list(x) for x in vector]
                dump_list_as_file(vector_average,args.react,args.single)
            else:
                print 'Transcript {} has insufficient values per bins'.format(args.single)
        except KeyError:
            print 'Your query key {} was not found in {}.'.format(args.single,args.react)
    
    #All Mode
    elif args.all_transcripts and not args.single and not args.multi:
        passed,failed,gummy = 0,0,[0]*args.bins
        for v in transcripts.values():
            vector = split_values(v,args.bins)
            if assert_bin_sizes(vector,args.size):
                vector_average = [average_sub_list(x) for x in vector]
                gummy = add_bins(gummy,vector_average)
                passed+=1
            else:
                failed+=1
        final_bin = [g/passed for g in gummy]
        dump_list_as_file(final_bin,args.react,'all')
        print '{} transcripts failed to have at least {} values in every bin of {} bins'.format(str(failed),str(args.size),str(args.bins))
    
    #Specified Mode
    elif args.multi and not args.single and not args.all_transcripts:
        targets = read_in_target_transcripts(args.multi)
        gummy,passed = [0]*args.bins,0
        for transcript in targets.keys():
            try:
                vector = split_values(transcripts[transcript],args.bins)
                if assert_bin_sizes(vector,args.size):
                    vector_average = [average_sub_list(x) for x in vector]
                    gummy = add_bins(gummy,vector_average)
                    passed+=1
                else:
                    print 'Transcript {} has insufficient values per bins, ignoring'.format(transcript)
            except KeyError:
                print 'Your query key {} was not found in {}.'.format(transcript,args.react)
        final_bin = [g/passed for g in gummy]
        dump_list_as_file(final_bin,args.react,args.multi.replace('.txt',''))

    #In case of Derp
    else:
        print 'Error. Arguments -single, -multi and -all_transcripts are mutually exclusive.'

if __name__ == '__main__':
    main()
