#!/usr/bin/env python

'''
Structure-Fold2 - react_windows.py
Compares two <.react> files, and outputs both a <.csv> and a <.fasta> of certain reactivity windows, i.e. gains, losses, or absolute differences.
The <.fasta> and <.react> files MUST be perfectly parallel in terms of sequence length. The <.fasta> can be put into meme, while the <.csv> may be
used for any purpose. It may be a good idea to run this program several times, i.e. once without any filters, then once again to get losses/gains/differences
and investigate all of these separately. 
'''

#Imports
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse

#Functions
def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

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

def subtract_reacts(reacts_A,reacts_B):
    '''subtracts B from A, i.e. A-B'''
    result = []
    spinner = zip(reacts_A,reacts_B)
    for bit in spinner:
        try:
            diff = bit[0]-bit[1]
            result.append(diff)
        except TypeError:
            result.append('NA')
    return result

def diff_reacts(reacts_A,reacts_B):
    '''Absolute Value of Difference, two lists'''
    result = []
    spinner = zip(reacts_A,reacts_B)
    for bit in spinner:
        try:
            diff = abs(bit[0]-bit[1])
            result.append(diff)
        except TypeError:
            result.append('NA')
    return result

def get_sum(alist):
    '''returns the sum of a list, removing NAs'''
    return sum([x for x in alist if x != 'NA'])

def gini(list_of_values):
    '''Returns the Gini value of a list of values.'''
    try:
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area
    except ZeroDivisionError:
        return 'NA'

def hot_stepper(fasta_seqs,control_reacts,experimental_reacts,window=50,step=20):
    '''Makes lists of reactivity windows'''
    disaster_area = []
    for k, v in fasta_seqs.items():
        try:
            control,experimental = control_reacts[k],experimental_reacts[k]
            change = subtract_reacts(experimental,control)
            abs_change = diff_reacts(experimental,control)
            value_windows = [change[i:i+window] for i in xrange(0, len(change)-(window-1), step)]
            abs_value_windows = [abs_change[i:i+window] for i in xrange(0, len(change)-(window-1), step)]
            seq_windows = [v[i:i+window] for i in xrange(0, len(v)-(window-1), step)]
            control_windows = [control[i:i+window] for i in xrange(0, len(v)-(window-1), step)]
            exp_windows = [experimental[i:i+window] for i in xrange(0, len(v)-(window-1), step)]
            if len(value_windows) > 1:
                for j in range(0, len(value_windows)):
                    sub_control,sub_exp = gini([z for z in control_windows[j] if z != 'NA']),gini([q for q in exp_windows[j] if q != 'NA'])
                    delta_gini = sub_exp - sub_control if all([ga != 'NA' for ga in [sub_control,sub_exp]]) else 'NA'
                    tidbit = tuple([get_sum(value_windows[j]),get_sum(abs_value_windows[j]),delta_gini,seq_windows[j],k,str(j),control_windows[j],exp_windows[j]])
                    disaster_area.append(tidbit)
        except KeyError:
            continue
    return disaster_area
    
def ultra_dumper(catastrophe,out='bits.csv'):
    '''Writes out a simple <.csv> file'''
    with open(out,'w') as g:
        g.write('net_change,total_change,gini_change,seq,transcript,step\n')
        for item in catastrophe:
            g.write(','.join([str(x) for x in item[0:6]])+'\n')

def write_out_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for item in info:
            name = '_'.join([item[4],item[5]])
            seq = item[3]
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n') 

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def fork_react_dicts(giant_info_list):
    '''Creates two sub-dictionaries in pseudo <.react> format'''
    x_control,x_experimetal = {},{}
    for item in giant_info_list:
        name = '_'.join([item[3],item[4]])
        control_frame, exp_frame = item[5],item[6]
        x_control[name] = control_frame
        x_experimetal[name] = exp_frame
    return x_control,x_experimetal

def write_react(react_info,outfyle='out.react'):
    '''Writes the <.react> back out'''
    with open(outfyle,'w') as g:
        for transcript, data in react_info.items():
            g.write(transcript+'\n')
            g.write('\t'.join([str(q) for q in data])+'\n')

def make_list_strict(mega_list):
    '''derp'''
    common_keys,new_list = {},[]
    while len(mega_list) !=0:
        item = mega_list[0]
        if item[3] in common_keys:
            del mega_list[0]
        else:
            common_keys[item[3]] = None
            new_list.append(item)
            del mega_list[0]
    return new_list


def write_react_fork(info,control_out,experimental_out):
    '''Writes two new reacts'''
    with open(control_out,'w') as g, open(experimental_out,'w') as h:
        for data in info:
            seq_name = '_'.join([data[3],data[4]])
            control_data = data[5]
            exp_data = data[6]
            g.write(seq_name+'\n')
            h.write(seq_name+'\n')
            g.write('\t'.join([str(q) for q in control_data])+'\n')
            h.write('\t'.join([str(z) for z in exp_data])+'\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Generate Windows of changing reactivity, will subtract control from experimental')
    parser.add_argument('control',type=str,help='control <.react> file')
    parser.add_argument('experimental',type=str,help='experimental <.react> file')
    parser.add_argument('fasta',type=str,help='<.fasta> to pull sequences from')
    parser.add_argument('-wlen',type=int, default=50, help='[default = 50] Window Length (bp)')
    parser.add_argument('-wstep',type=int, default=20, help='[default = 20] Window Step (bp)')
    parser.add_argument('-outname',type=str,default=None, help='Change the name of the outfile, overrides default')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-sort',type=str.upper, default = None,choices = ['L','G','D'],help = 'Sort windows by Loss/Gain/Delta')
    parser.add_argument('-perc',type=int, default=100, help='[default = 100] Percentage of sorted windows to retain')
    parser.add_argument('-singular',action='store_true',default=False,help = 'Limit to top window per transcript after sorting')
    parser.add_argument('-strict',action='store_true',default=False,help = '-sort removes non gain/loss/delta before perc filtering')
    parser.add_argument('-fastaout',action='store_true',default=False,help = 'Write out windows in <.fasta> format as well')
    parser.add_argument('-reactout',action='store_true',default=False,help = 'Write out windows in <.react> format as well')
    args = parser.parse_args()
    
    #Generate Outfile name
    default_name = '_'.join([args.control.replace('.react',''),args.experimental.replace('.react',''),str(args.wlen)+'win',str(args.wstep)+'step'])+'.csv'
    out_name = default_name if args.outname == None else check_extension(args.outname,'.csv')
    
    #Read in both groups of reactivities, fasta file with sequences
    control_reactivty,experimental_reactivty = read_reactivities(args.control),read_reactivities(args.experimental)
    target_seqs = read_fasta(args.fasta)
    
    #Truncate seqs to those with good coverage in both conditions if user provides a list
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(target_seqs,covered)
    
    #Generate Windows
    loud_noises = hot_stepper(target_seqs,control_reactivty,experimental_reactivty,args.wlen,args.wstep)

    #Apply filter
    if args.sort == 'L':
        loud_noises = loud_noises if args.strict == False else [item for item in loud_noises if item[0] < 0]
        loud_noises = sorted(loud_noises,reverse=False)[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'loss.csv')

    if args.sort == 'G':
        loud_noises = loud_noises if args.strict == False else [item for item in loud_noises if item[0] > 0]
        loud_noises = sorted(loud_noises,reverse=True)[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'gain.csv')
    
    if args.sort == 'D':
        loud_noises = loud_noises if args.strict == False else [item for item in loud_noises if item[0] != 0]
        loud_noises = sorted(loud_noises,reverse=True,key=lambda x:x[1])[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'delta.csv')
    
    #Singular filter if applicable
    if args.sort and args.singular:
        loud_noises = make_list_strict(loud_noises)

    #Output Suite
    #Output for <.csv>
    ultra_dumper(loud_noises,out_name)
    
    #Output for <.fasta>
    if args.fastaout == True:
        write_out_fasta(loud_noises,out_name.replace('.csv','.fasta'))
    
    #Output for <.react>
    if args.reactout == True:
        new_control_file = '_'.join([args.control.replace('.react',''),str(args.wlen)+'win',str(args.wstep)+'step'])+'.react'
        new_exp_file = '_'.join([args.experimental.replace('.react',''),str(args.wlen)+'win',str(args.wstep)+'step'])+'.react'
        if args.sort == 'L':
            new_control_file = new_control_file.replace('.react','_'+str(args.perc)+'loss.react')
            new_exp_file = new_exp_file.replace('.react','_'+str(args.perc)+'loss.react')
        if args.sort == 'G':
            new_control_file = new_control_file.replace('.react','_'+str(args.perc)+'gain.react')
            new_exp_file = new_exp_file.replace('.react','_'+str(args.perc)+'gain.react')
        if args.sort == 'D':
            new_control_file = new_control_file.replace('.react','_'+str(args.perc)+'delta.react')
            new_exp_file = new_exp_file.replace('.react','_'+str(args.perc)+'delta.react')
        write_react_fork(loud_noises,new_control_file,new_exp_file)

if __name__ == '__main__': 
    main()
