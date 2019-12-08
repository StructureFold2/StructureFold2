#!/usr/bin/env python2

'''
This script takes two <.rtsc> along with the base <.fasta> to calculate the reactivity scores of each nucelotide.
If not given a <.scale> file, the script will generate a <.scale> to be used with this script when calculating any reactvities you wish to do a comparison to.
Transcripts which produce a 0 average on the normalization scale will not have their final reactivity calculated, and these transcripts may be logged to a file (optional).
If given a <.scale> file, it will apply it as the normalization scale. Transcripts not included in the <.scale> will not be calculated and can be logged to a file (optional).
Normalization (2-8%) may be turned off for the reactivity calculation, as can taking the natural log of the reactivity values, via options.
'''

#Imports
import math
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from itertools import islice

#Functions
def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary, transcript_name:transcript_sequence'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary, transcript_name:[list of stop numbers]'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) for x in stops.split('\t')]
    return information

def calculate_raw_reactivity(reagent_minus,reagent_plus,nlog_off=False):
    '''Calculates raw reactivity, with or without the natural log'''
    key_set,data_out = set(reagent_plus.keys()).intersection(set(reagent_minus.keys())),{}
    for key in key_set:
        plus_vector = [math.log(value+1,math.e) for value in reagent_plus[key]] if nlog_off == False else reagent_plus[key]
        minus_vector = [math.log(value+1,math.e) for value in reagent_minus[key]] if nlog_off == False else reagent_minus[key]
        sum_plus,sum_minus,length = sum(plus_vector),sum(minus_vector),len(plus_vector)
        if sum_plus != 0 and sum_minus != 0:
            nrm_plus_vector = [float(y)/float(sum_plus)*length for y in plus_vector]
            nrm_minus_vector = [float(z)/float(sum_minus)*length for z in minus_vector]
            plus_minus = [max(0,a-b) for a, b in zip(nrm_plus_vector,nrm_minus_vector)]
            data_out[key] = plus_minus
    return data_out

def generate_normalization_scale(derived_reactivities,transcript_seqs,specificity):
    '''Generates the 2-8% scale to normalize against'''
    data = {}
    for transcript, reactivities in derived_reactivities.items():
        sequence = transcript_seqs[transcript]
        accepted = sorted([reactivities[k] for k in range(1,len(reactivities)) if sequence[k-1] in specificity],reverse=True)
        top = accepted[int(len(accepted)*0.02):int(len(accepted)*0.1)]
        top_average = sum(top)/len(top) if len(top) > 0 else 0
        if top_average > 0:
            data[transcript] = top_average
    return data

def read_normalization_scale(normalization_file):
    '''Reads in a normalization scale file'''
    info = {}
    with open(normalization_file, 'r') as f:
        for line in f:
            if line.startswith('transcript'):
                continue
            else:
                transcript,value = line.strip().split(',')
                info[transcript] = float(value)
    return info

def write_normalization_scale(scale_dictionary,outfile):
    '''Writes out a normalization scale file'''
    with open(outfile,'w') as g:
        g.write(','.join(['transcript','value'])+'\n')
        for transcript, value in scale_dictionary.items():
            g.write(','.join([transcript,str(value)])+'\n')

def calculate_final_reactivity(derived_reactivities,transcript_sequences,specificity,threshold,nrm_scale,norm_off=False):
    '''Calculates the final reactivity'''
    data_out,missing_transcripts = {},{}
    for transcript, reactivities in derived_reactivities.items():
        if transcript in nrm_scale:
            normalizer = nrm_scale[transcript] if norm_off == False else 1
            sequence = transcript_sequences[transcript]
            normalized_values =[str(float('%.3f'%min((reactivities[x]/normalizer), threshold))) if sequence[x-1] in specificity else 'NA' for x in range(1,len(reactivities))]+['NA']
            data_out[transcript] = normalized_values
        else:
            missing_transcripts[transcript] = None
    return data_out,missing_transcripts

def write_out_reactivity_file(info, outfyle):
    '''Writes out a <.react> file'''
    with open(outfyle,'w') as g:
        for name, values in info.items():
            g.write(name+'\n')
            g.write('\t'.join(values)+'\n')

def write_out_missing_file(info,outfyle):
    '''Writes out the transcripts which were missing from the normalizaition scale.'''
    with open(outfyle,'w') as g:
        for transcript in info.keys():
            g.write(transcript+'\n')

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def main():
    parser = argparse.ArgumentParser(description='Generates <.react> files from two <.rtsc> files')
    parser.add_argument('control',type=str,help='Control <.rtsc> file')
    parser.add_argument('treatment',type=str,help='Reagent <.rtsc> file')
    parser.add_argument('fasta',type=str,help='Transcript <.fasta> file')
    parser.add_argument('-threshold',type=float,default=7.0,help='[default = 7.0] Reactivity Cap')
    parser.add_argument('-ln_off',action='store_true', help='Do not take the natural log of the stop counts')
    parser.add_argument('-nrm_off',action='store_true',help='Turn off 2-8'+u"\uFF05"+' normalization of the derived reactivity')
    parser.add_argument('-save_fails',action='store_true',help='Log transcripts with zero or missing scales')
    parser.add_argument('-scale',type=str,default=None, help='Provide a normalizaiton <.scale> for calculation')
    parser.add_argument('-bases',type=str,default='AC', help='[default = AC] Reaction Specificity, (AGCT) for SHAPE')
    parser.add_argument('-name',type=str,default=None, help='Change the name of the outfile, overrides default')
    parser.add_argument('-restrict',default = None, help = 'Limit analysis to these specific transcripts <.txt> ')
    args = parser.parse_args()
    
    #Create output name
    base_name = [x.split(os.sep)[-1].replace('.rtsc','') for x in [args.control,args.treatment]]
    log_tag = ['ln'] if args.ln_off==False else []
    nrm_tag = ['nrm'] if args.ln_off==False else []
    out_name = '_'.join(base_name+log_tag+nrm_tag)+'.react' if args.name == None else check_extension(args.name,'.react')

    #Read in data, apply restrictions if applicable
    control_data,treatment_data = read_in_rtsc(args.control),read_in_rtsc(args.treatment)
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(control_data,covered)
        filter_dictonary(treatment_data,covered)
    
    #Calculate Derived Reactivity
    data = calculate_raw_reactivity(control_data,treatment_data,args.ln_off)
    
    #Read in transcript sequences
    seqs = read_in_fasta(args.fasta)
    
    #Generate and write scale, or read a <.scale> file in
    normalizaiton_scale = generate_normalization_scale(data,seqs,args.bases) if args.scale == None else read_normalization_scale(args.scale)
    if args.scale == None:
        write_normalization_scale(normalizaiton_scale,out_name.replace('.react','.scale'))
    
    #Calculate Final Reactivity
    out_reactivity,out_missing = calculate_final_reactivity(data,seqs,args.bases,args.threshold,normalizaiton_scale,args.nrm_off)
    
    #Write Out
    write_out_reactivity_file(out_reactivity,out_name)
    if args.save_fails:
        write_out_missing_file(out_missing,out_name.replace('.react','_unresolvable_transcripts.txt'))

if __name__ == '__main__':
    main()
