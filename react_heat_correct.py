#!/usr/bin/env python

'''
Must be used on .react files that have already been filtered to only 
include the transcripts above the accepted coverage threshold, and must be the exact same transcripts!
'''

#Imports
import itertools
import argparse
import sys

#Functions
def read_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(itertools.islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!= 'NA' else 'NA' for x in reactivities.split()]
    return information

def write_react(react_info,outfyle='out.react'):
    '''Writes the <.react> back out'''
    with open(outfyle,'w') as g:
        for transcript, data in react_info.items():
            g.write(transcript+'\n')
            g.write('\t'.join([str(q) for q in data])+'\n')

def sum_react(react_dict):
    '''Sum of all the reactivities'''
    return sum([sum(filter(lambda x: isinstance(x, float), v)) for v in react_dict.values()])

def apply_correction(react_dict,correction):
    '''Applies a correction'''
    atarashii = {}
    for k, v in react_dict.items():
        atarashii[k] = [x*correction if x!= 'NA' else 'NA' for x in v]
    return atarashii

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects temperature differential derived reacts')
    parser.add_argument('lower',type=str,help='lower temp <.react> file')
    parser.add_argument('higher',type=str,help='higher temp <.react> file')
    parser.add_argument('-suffix',type=str,default='corrected',help='[default = corrected] Suffix for out files')
    args = parser.parse_args()
    
    #Workflow
    
    #Sum all reactivities
    cold_react,hot_react = map(read_reactivities,[args.lower,args.higher])
    cold_sum,hot_sum = map(sum_react,[cold_react,hot_react])
    
    if not cold_react.keys() == hot_react.keys():
        print 'Warning! Non-parallel transcript sets between reacts! Quitting.. '
        sys.exit() 

    #Calculate Corrections
    heat_correction = (hot_sum+cold_sum)/(2*hot_sum)
    cold_correction = (hot_sum+cold_sum)/(2*cold_sum)
    
    print 'Higher Temp values to be downscaled by factor: {}'.format(heat_correction)
    print 'Lower Temp values to be upscaled by factor: {}'.format(cold_correction)
    
    #Apply corrections
    new_hot = apply_correction(hot_react,heat_correction)
    new_cold = apply_correction(cold_react,cold_correction)
    
    #Write Out
    write_react(new_cold,args.lower.replace('.react','_'+args.suffix+'.react'))
    write_react(new_hot,args.higher.replace('.react','_'+args.suffix+'.react'))
    

if __name__ == '__main__': 
    main()
