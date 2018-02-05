#!/usr/bin/env python

'''
Info soon!
'''

#Imports
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import re

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

def generate_motif_coords(sequence,motif):
    '''Returns all coordinates of matching patterns'''
    return [(match.start(), match.end()) for match in re.finditer(motif, sequence)]

def base_pull(item,index,fill='-'):
    '''Pulls from phantom areas'''
    if index < 0:
        thing = fill
    else:
        try:
            thing = item[index]
        except IndexError:
            thing = fill
    return thing

def buffered_pull(target,coord_tupe,fpbuffer=0,tpbuffer=0,e_fill='-'):
    '''Pulls the regions, fills in either side if requested window is out of range with base_pull'''
    tp_side = [base_pull(target,x,e_fill) for x in range(coord_tupe[1],coord_tupe[1]+tpbuffer)]
    fp_side = [base_pull(target,x,e_fill) for x in range(coord_tupe[0]-fpbuffer,coord_tupe[0])]
    if type(target) == str:
        return fp_side+list(target)[coord_tupe[0]:coord_tupe[1]]+tp_side
    elif type(target) == list:
        return fp_side+target[coord_tupe[0]:coord_tupe[1]]+tp_side
    else:
        return None

def cold_stepper(fasta_seqs,control_reacts,experimental_reacts,motif,fpbuffer,tpbuffer):
    '''Does the thing'''
    sanity_area = {}
    for k, v in fasta_seqs.items():
        try:
            control,experimental = control_reacts[k],experimental_reacts[k]
            change = subtract_reacts(experimental,control)
            indexes = generate_motif_coords(v,motif)
            for i in indexes:
                new_seq = buffered_pull(v,i,fpbuffer,tpbuffer,'-')
                con_numbs = buffered_pull(control,i,fpbuffer,tpbuffer,'NA')
                exp_numbs = buffered_pull(experimental,i,fpbuffer,tpbuffer,'NA')
                chg_nubms = buffered_pull(change,i,fpbuffer,tpbuffer,'NA')
                #q_key = '_'.join([k, '~'.join([str(i[0]+1),str(i[1]+1)])])
                q_key = (k,str(i[0]+1),str(i[1]+1))
                sanity_area[q_key] = [new_seq,con_numbs,exp_numbs,chg_nubms]
        except KeyError:
            continue
    return sanity_area

def dump_csv(information,outfile,motif,fpbuffer,tpbuffer):
    '''Writes out the informaiton to a <.csv>'''
    #generate header
    base_header = ['transcript','motif_start','motif_end']
    z_length = fpbuffer+len(motif)+tpbuffer
    control_columns = ['_'.join(['c',str(x)]) for x in range(1,z_length+1)]
    experimental_columns = ['_'.join(['e',str(x)]) for x in range(1,z_length+1)]
    delta_columns = ['_'.join(['d',str(x)]) for x in range(1,z_length+1)]
    FP_columns = ['_'.join(['fp',str(x)]) for x in range(1,fpbuffer+1)]
    M_colmuns = ['_'.join(['m',str(x)]) for x in range(1,len(motif)+1)]
    TP_columns = ['_'.join(['tp',str(x)]) for x in range(1,tpbuffer+1)]
    seq_columns = FP_columns+M_colmuns+TP_columns
    EX_columns = ['motif_Change','motif_Increase','motif_Decrease','motif_Delta']
    Full_Header = base_header+seq_columns+control_columns+experimental_columns+delta_columns+EX_columns
    #
    with open(outfile,'w') as g:
        g.write(','.join(Full_Header)+'\n')
        for k, v in sorted(information.items()):
            base = list(k)
            for item in v:
                base.extend(item)
            #
            z_numbers = [x for x in v[3][fpbuffer:len(v)-tpbuffer] if x != 'NA']
            change = sum(z_numbers)
            positive = sum([n for n in z_numbers if n > 0])
            negative = sum([n for n in z_numbers if n < 0])
            delta = sum([abs(n) for n in z_numbers])
            EX_Vector = [change,positive,negative,delta]
            base.extend(EX_Vector)
            #
            base = [str(z) for z in base]
            g.write(','.join(base)+'\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Searches and returns reactivity vectors for a targeted motif')
    parser.add_argument('control',type=str,help='control <.react> file')
    parser.add_argument('experimental',type=str,help='experimental <.react> file')
    parser.add_argument('fasta',type=str,help='<.fasta> to pull sequences from')
    parser.add_argument('motif',type=str,help='nucleotide motif')
    parser.add_argument('-fp',default=5,type=int, help='Bases to include 5\' of the motif')
    parser.add_argument('-tp',default=5,type=int, help='Bases to include 3\' of the motif')
    parser.add_argument('-outname',type=str,default=None, help='Change the name of the outfile, overrides default')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    args = parser.parse_args()
    
    #Generate Name
    default_name = '_'.join([args.control.replace('.react',''),args.experimental.replace('.react',''),args.motif,str(args.fp)+'fp',str(args.tp)+'tp'])+'.csv'
    out_name = default_name if args.outname == None else args.outname
    
    #Read in both groups of reactivities, fasta file with sequences
    control_reactivty,experimental_reactivty = read_reactivities(args.control),read_reactivities(args.experimental)
    target_seqs = read_fasta(args.fasta)
    
    #Truncate seqs to those with good coverage in both conditions if user provides a list
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(target_seqs,covered)

    #Pop windows
    quiet_noises = cold_stepper(target_seqs,control_reactivty,experimental_reactivty,args.motif,args.fp,args.tp)
    
    #Write
    dump_csv(quiet_noises,out_name,args.motif,args.fp,args.tp)



if __name__ == '__main__': 
    main()
