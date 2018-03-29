#!/usr/bin/env python

'''
Details Net/Negative/Positive change along a set number of bases from either the 5' or 3' end, from one or more
transcripts between two <.react> files.
'''

#Imports
import argparse
import os
from itertools import islice
from collections import Counter
import numpy

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

def sd_counter(acounter):
    ''''''
    new_counter = Counter()
    for index,sub_counter in acounter.items():
        new_counter[index] = numpy.std([derp for herp in [[k]*v for k, v in sub_counter.items()] for derp in herp])
    return new_counter

def reactivity_pattern_delta(reacts_a,reacts_b,positions=150,t_prime=False,raw=False,sd=False):
    '''Generates distibutions of change between two sets of reactivity values'''
    delta,positive,negative, = Counter(),Counter(),Counter()
    delta_norm,pos_norm,neg_norm = Counter(),Counter(),Counter()
    change,change_norm = Counter(),Counter()
    c_prior,e_prior = Counter(),Counter()
    c_norm,e_norm = Counter(),Counter()
    delta_sd,change_sd,pos_sd,neg_sd,c_sd,e_sd = {},{},{},{},{},{}
    
    #
    transcripts = set(reacts_a.keys()).intersection(set(reacts_b.keys()))
    for transcript in transcripts:
        a,b = reacts_a[transcript],reacts_b[transcript]
        try:
            values = [(i+1,b[i]-a[i]) for i in range(0,len(a))][:positions] if not t_prime else [(i-len(a),b[i]-a[i]) for i in range(0,len(a))][-positions:]
            c_values = [(i+1,a[i]) for i in range(0,len(a))][:positions] if not t_prime else [(i-len(a),a[i]) for i in range(0,len(a))][-positions:]
            e_values = [(i+1,b[i]) for i in range(0,len(b))][:positions] if not t_prime else [(i-len(b),b[i]) for i in range(0,len(b))][-positions:]
            #
            for item in values:
                delta[item[0]]+= item[1]
                delta_norm[item[0]]+=1
                delta_sd.setdefault(item[0], Counter())[item[1]]+=1
                if item[1] !=0:
                    change[item[0]]+= item[1]
                    change_norm[item[0]]+=1
                    change_sd.setdefault(item[0], Counter())[item[1]]+=1
                if item[1] > 0:
                    positive[item[0]]+= item[1]
                    pos_norm[item[0]]+=1
                    pos_sd.setdefault(item[0], Counter())[item[1]]+=1
                if item[1] < 0:
                    negative[item[0]]+= item[1]
                    neg_norm[item[0]]+=1
                    neg_sd.setdefault(item[0], Counter())[item[1]]+=1
            #
            for c_item in c_values:
                c_prior[c_item[0]]+=c_item[1]
                c_norm[c_item[0]]+=1
                c_sd.setdefault(c_item[0], Counter())[c_item[1]]+=1
            for e_item in e_values:
                e_prior[e_item[0]]+=e_item[1]
                e_norm[e_item[0]]+=1
                e_sd.setdefault(e_item[0], Counter())[e_item[1]]+=1
        except IndexError:
            continue
    #Permute to type.
    if raw == False and sd==False:
        delta,positive,negative = average_counter(delta,delta_norm),average_counter(positive,pos_norm),average_counter(negative,neg_norm)
        change,control_base,experimental_base = average_counter(change,change_norm),average_counter(c_prior,c_norm),average_counter(e_prior,e_norm)
    if raw == False and sd==True:
        delta,positive,negative = sd_counter(delta_sd),sd_counter(pos_sd),sd_counter(neg_sd)
        change,control_base,experimental_base = sd_counter(change_sd),sd_counter(c_sd),sd_counter(e_sd)
    if raw == True:
        control_base,experimental_base = c_prior,e_prior
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
        if coordinate not in control_base:
            control_base[coordinate] = 0.0
        if coordinate not in experimental_base:
            experimental_base[coordinate] = 0.0
    #Return items
    return (delta,change,positive,negative,control_base,experimental_base)

def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = None
    return good

def write_out_results(data_tuple,use_raw,outfyle='out.csv'):
    '''Writes out the data in a <.csv>'''
    base_header = ['Position','Control','Experimental','Average_Delta','Average_Change','Average_Postive_Change','Average_Negative_Change']
    base_header = base_header if use_raw == False else [x.replace('Average','Raw') for x in base_header]
    with open(outfyle,'w') as g:
        g.write(','.join(base_header)+'\n')
        delta,change,pos,neg,cbase,ebase = data_tuple
        for key in sorted(delta.keys()):
            x_line = ','.join([str(key),str(cbase[key]),str(ebase[key]),str(delta[key]),str(change[key]),str(pos[key]),str(neg[key])])+'\n'
            g.write(x_line)

def write_out_results_as_bins(data_tuple,use_raw,outfyle='out.csv',nbin=10):
    '''De-Derp the thingies'''
    base_header = ['Position','Control','Experimental','Average_Delta','Average_Change','Average_Postive_Change','Average_Negative_Change']
    base_header = base_header if use_raw == False else [x.replace('Average','Raw') for x in base_header]
    with open(outfyle,'w') as g:
        g.write(','.join(base_header)+'\n')
        delta,change,pos,neg,cbase,ebase = data_tuple
        sub_bins = [sorted(delta.keys())[x:x+nbin] for x in range(0,len(delta.keys()),nbin)]
        new_keys = ['-'.join([str(x) for x in [sublist[0],sublist[-1]]]) for sublist in sub_bins ]
        for key, box in zip(new_keys,sub_bins):
            data_vector = [sum([goober[index] for index in box])/len(box) for goober in [cbase,ebase,delta,change,pos,neg]]
            z_line = ','.join([str(key)]+[str(v) for v in data_vector])+'\n'
            g.write(z_line)

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
    parser.add_argument('-sd',action='store_true', help='Report standard deviation rather than average')
    parser.add_argument('-raw',action='store_true', help='Report raw sums, overrides -sd')
    parser.add_argument('-single', type=str,default=None, help='Target transcript to bin from the file')
    parser.add_argument('-multi', type=str,default=None, help='Flat <.txt> of transcripts to process, one per line')
    parser.add_argument('-all_transcripts',action='store_true', help='Use all transcripts in the <.react> file')
    parser.add_argument('-bin_numb',type=int, default=0, help='Bin all results by input')
    args = parser.parse_args()
    
    #Read in the <.react> files
    control,experimental = read_in_derived_reactivities(args.control),read_in_derived_reactivities(args.experimental)
    name_bit = 'FP' if args.tp == False else 'TP'
    type_bit = 'Average' if args.sd == False else 'SD'
    type_bit = type_bit if args.raw == False else 'Raw'

    #Single Mode
    if args.single and not args.all_transcripts and not args.multi:
        try:
            c_tiny,e_tiny = {},{}
            c_tiny[args.single] = control[args.single]
            e_tiny[args.single] = experimental[args.single]
            out_data = reactivity_pattern_delta(c_tiny,e_tiny,args.n,args.tp,args.raw,args.sd)
            base_name = [x.split(os.sep)[-1].replace('.react','') for x in [args.control,args.experimental]]
            if args.bin_numb == 0:
                outfyle = '_'.join(base_name+[args.single,name_bit,type_bit,str(args.n)+'bp'])+'_delta.csv'
                write_out_results(out_data,args.raw,outfyle)
            else:
                outfyle = '_'.join(base_name+[args.single,name_bit,type_bit,str(args.n)+'bp',str(args.bin_numb)+'bin'])+'_delta.csv'
                write_out_results_as_bins(out_data,args.raw,outfyle,args.bin_numb)
        except KeyError:
            print 'Your query key {} was not found in {}.'.format(args.single,args.react)
    
    #All Transcripts Mode
    elif args.all_transcripts and not args.single and not args.multi:
        out_data = reactivity_pattern_delta(control,experimental,args.n,args.tp,args.raw,args.sd)
        base_name = [x.split(os.sep)[-1].replace('.react','') for x in [args.control,args.experimental]]
        if args.bin_numb == 0:
            outfyle = '_'.join(base_name+['all',name_bit,type_bit,str(args.n)+'bp'])+'_delta.csv'
            write_out_results(out_data,args.raw,outfyle)
        else:
            outfyle = '_'.join(base_name+['all',name_bit,type_bit,str(args.n)+'bp',str(args.bin_numb)+'bin'])+'_delta.csv'
            write_out_results_as_bins(out_data,args.raw,outfyle,args.bin_numb)

    #Specified Mode
    elif args.multi and not args.single and not args.all_transcripts:
        keepers = read_in_target_transcripts(args.multi)
        clean_dictionary(keepers,control)
        clean_dictionary(keepers,experimental)
        out_data = reactivity_pattern_delta(control,experimental,args.n,args.tp,args.raw,args.sd)
        base_name = [x.split(os.sep)[-1].replace('.react','') for x in [args.control,args.experimental]]
        if args.bin_numb == 0:
            outfyle = '_'.join(base_name+[args.multi.split('.')[0],name_bit,type_bit,str(args.n)])+'_delta.csv'
            write_out_results(out_data,args.raw,outfyle)
        else:
            outfyle = '_'.join(base_name+[args.multi.split('.')[0],name_bit,type_bit,str(args.n)+'bp',str(args.bin_numb)+'bin'])+'_delta.csv'
            write_out_results_as_bins(out_data,args.raw,outfyle,args.bin_numb)
    
    #In case of Derp
    else:
        print 'Error. Arguments -single, -multi and -all_transcripts are mutually exclusive.'

if __name__ == '__main__':
    main()
