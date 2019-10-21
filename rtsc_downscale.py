#!/usr/bin/env python2

'''
Downscales <.rtsc> files, keeping a given percentage of the original RT hits.
This downscaling cab be done statically, simply multiplying every RT count by a specified percentage
This downscaling can be done pseudo-randomly, removing hits randomly down to a specified percentage.

#ISSUE - because of the way we are randomly removing, bases with lots of stops aren't losing anything
and bases with few are being wiped entirely. Should change to an implementation where at each base, count
the hits and each hit on each base has X percent change of retention, not throwing X total craters at it.
'''

#Imports
import itertools
import random
import glob
import argparse

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

def keep_static_percentage(alist,percentage=.50):
    '''Evenly retains a given percentage of RT hits. Can generate fractional hits.'''
    return [float(i)*percentage for i in alist]

def keep_random_downsample(num_list,percentage=.50):
    '''Randomly retains a given percentage of RT hits'''
    total,length,internal = sum(num_list),len(num_list),num_list[:]
    keeps = int(total*percentage)
    removes = total-keeps
    while removes > 0:
        target = random.randrange(length)
        if internal[target] > 0:
            internal[target]-=1
            removes+=-1
    return internal

def random_per_stop(alist,chance=.50):
    '''Roll a Die for each RT stop at each positon'''
    newlist,new_value = [],int(chance*1000)
    for value in alist:
        dice_rolls = [random.randint(1,1000) for stop in range(0,value)]
        keeps = filter(lambda x:x <= new_value, dice_rolls)
        newlist.append(len(keeps))
    return newlist

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

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def main():
    parser = argparse.ArgumentParser(description='Determines approximate transcript abundance based on <.rtsc> files.')
    parser.add_argument('-f',default=None, nargs='+', help = '<.rtsc> to operate on')
    parser.add_argument('mode',type=str.upper, choices = ['PERCENT','RANDOMREAD','RANDOMPOSITION'])
    parser.add_argument('-ratio',type=float,default=.50, help='[default = 0.50] Fraction of RT stops to retain')
    parser.add_argument('-restrict',default = None, help = 'Limit analysis to these specific transcripts <.txt> ')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    args = parser.parse_args()

    #Files to operate on, dictionary of functions
    fyle_lyst = sorted(glob.glob('*.rtsc')) if args.f == None else sorted(args.f)
    downscale_methods = {'PERCENT':keep_static_percentage,'RANDOMPOSITION':keep_random_downsample,'RANDOMREAD':random_per_stop}

    #Iterate through file(s) 
    for fyle in fyle_lyst:
        data,new_fyle = read_in_rtsc(fyle),fyle.replace('.rtsc','_'+args.mode+'_'+str(args.ratio).replace('.','')+'.rtsc')
        
        #You could be doing random for a long time if you do not filter by coverage.
        if args.restrict != None:
            covered = get_covered_transcripts(args.restrict)
            data = dict([(x,z) for x,z in data.items() if x in covered])

        new_data = dict([(k, downscale_methods[args.mode](v,args.ratio)) for k, v in data.items()])
        write_rtsc(new_data,args.sort,new_fyle)

#Main
if __name__ == '__main__':
    main()
