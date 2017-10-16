#!/usr/bin/python

'''
Details Net/Negative/Positive change along a set number of bases from either the 5' or 3' end, from one or more
transcripts between two <.react> files.
'''

#Imports
import argparse
from itertools import islice
from collections import Counter

#Functions
def read_in_derived_reactivities(reactvity_fyle):
    '''Reads in a reactivity file'''
    information = {}
    with open(reactvity_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!='NA' else 0.0 for x in reactivities.split('\t') ]
    return information

def average_counter(counter_1,counter_2):
    '''Divides all values of 1 by values in 2'''
    new_counter = Counter()
    for key, value in counter_1.items():
        new_counter[key] = value/counter_2[key]
    return new_counter

def reactivity_pattern_delta(reacts_a,reacts_b,positions=150,t_prime=False,raw=False):
    '''Generates distibutions of change between two sets of reactivity values'''
    delta,positive,negative, = Counter(),Counter(),Counter()
    delta_norm,pos_norm,neg_norm = Counter(),Counter(),Counter()
    change,change_norm = Counter(),Counter()
    #
    transcripts = set(reacts_a.keys()).intersection(set(reacts_b.keys()))
    for transcript in transcripts:
        a,b = reacts_a[transcript],reacts_b[transcript]
        try:
            values = [(i+1,b[i]-a[i]) for i in range(0,len(a))][:positions] if not t_prime else [(i-len(a),b[i]-a[i]) for i in range(0,len(a))][-positions:]
            for item in values:
                delta[item[0]]+= item[1]
                delta_norm[item[0]]+=1
                if item[1] !=0:
                    change[item[0]]+= item[1]
                    change_norm[item[0]]+=1
                if item[1] > 0:
                    positive[item[0]]+= item[1]
                    pos_norm[item[0]]+=1
                if item[1] < 0:
                    negative[item[0]]+= item[1]
                    neg_norm[item[0]]+=1
        except IndexError:
            continue
    #Average if enabled.
    if not raw:
        delta,positive,negative = average_counter(delta,delta_norm),average_counter(positive,pos_norm),average_counter(negative,neg_norm)
        change = average_counter(change,change_norm)
    #Check for gaps, then fill in zeros
    for coordinate in delta_norm.keys():
        if coordinate not in delta:
            delta[coordinate] = 0.0
        if coordinate not in positive:
            positive[coordinate] = 0.0
        if coordinate not in negative:
            negative[coordinate] = 0.0
        if coordinate not in change:
            change[coordinate] = 0.0
    #Return items
    return (delta,change,positive,negative)

def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = None
    return good

def write_out_results(data_tuple,outfyle='out.csv'):
    '''Writes out the data in a <.csv>'''
    with open(outfyle,'w') as g:
        g.write(','.join(['Position','Average_Delta','Average_Change','Average_Positive_Change','Average_Negative_Change'])+'\n')
        delta,change,pos,neg = data_tuple
        for key in sorted(delta.keys()):
            x_line = ','.join([str(key),str(delta[key]),str(change[key]),str(pos[key]),str(neg[key])])+'\n'
            g.write(x_line)

def clean_dictionary(keep_values,target_dict):
    '''removes values'''
    for k in target_dict.keys():
        if k not in keep_values:
            del target_dict[k]

def main():
    parser = argparse.ArgumentParser(description='Generates changes along n bp from either end of one or more transcripts')
    parser.add_argument('control', type=str, help='Control <.react> file')
    parser.add_argument('experimental',type=str, help='Experimental <.react> file')
    parser.add_argument('-n',type=int,default=100, help='Number of bases, default 100')
    parser.add_argument('-tp',action='store_true', help='Start from the 3\' end (5\' default)')
    parser.add_argument('-raw',action='store_true', help='Do not average changes')
    parser.add_argument('-single', type=str,default=None, help='Target transcript to bin from the file')
    parser.add_argument('-multi', type=str,default=None, help='Flat <.txt> of transcripts to bin, one per line')
    parser.add_argument('-all_transcripts',action='store_true', help='Use all transcripts in the <.react> file')
    args = parser.parse_args()
    
    #Read in the <.react> files
    control,experimental = read_in_derived_reactivities(args.control),read_in_derived_reactivities(args.experimental)
    name_bit = 'FP' if args.tp == False else 'TP'
    
    #Single Mode
    if args.single and not args.all_transcripts and not args.multi:
        try:
            c_tiny,e_tiny = {},{}
            c_tiny[args.single] = control[args.single]
            e_tiny[args.single] = experimental[args.single]
            out_data = reactivity_pattern_delta(c_tiny,e_tiny,args.n,args.tp,args.raw)
            outfyle = '_'.join([args.control.replace('.react',''), args.experimental.replace('.react',''),args.single,name_bit])+'_delta.csv'
            write_out_results(out_data,outfyle)
        except KeyError:
            print 'Your query key {} was not found in {}.'.format(args.single,args.react)
    
    #All Mode
    elif args.all_transcripts and not args.single and not args.multi:
        out_data = reactivity_pattern_delta(control,experimental,args.n,args.tp,args.raw)
        outfyle = '_'.join([args.control.replace('.react',''), args.experimental.replace('.react',''),'all',name_bit])+'_delta.csv'
        write_out_results(out_data,outfyle)

    #Specified Mode
    elif args.multi and not args.single and not args.all_transcripts:
        keepers = read_in_target_transcripts(args.multi)
        clean_dictionary(keepers,control)
        clean_dictionary(keepers,experimental)
        out_data = reactivity_pattern_delta(control,experimental,args.n,args.tp,args.raw)
        outfyle = '_'.join([args.control.replace('.react',''), args.experimental.replace('.react',''),args.multi.split('.')[0],name_bit])+'_delta.csv'
        write_out_results(out_data,outfyle)

    #In case of Derp
    else:
        print 'Error. Arguments -single, -multi and -all_transcripts are mutually exclusive.'

if __name__ == '__main__':
    main()
