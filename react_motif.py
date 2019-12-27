#!/usr/bin/env python2

#Imports
import os
import argparse
import itertools
from sf2libs import structure_io
from sf2libs.react_utils import MotifReport

def read_motif(user_input):
    '''Reads either a file or a motif, returns list of motifs'''
    if os.access(user_input,os.R_OK):
        with open(user_input,'r') as f:
            basic_motifs = [line.strip().upper() for line in f]
            return list(set(basic_motifs))
    else:
        return [user_input.upper()]
 
def write_motif_csv(motif_report,outfile='file.csv'):
    '''Writes out the default output'''
    row_names = ['transcript','query','start','end']
    m_names = ['m'+str(i) for i in range(1,len(motif_report.motif)+1)]
    fp_names = ['fp'+str(q) for q in range(1,motif_report.FP+1)]
    tp_names = ['tp'+str(p) for p in range(1,motif_report.TP+1)]
    val_types = ['A','B','D']
    value_names = [lab+str(n) for lab in val_types for n in range(1,motif_report.size+1)]
    prefixes,stats_types = ['query','fp','tp'],['change','increases','decreases','delta']
    stats_names = [lemon+'_'+lime for lemon in prefixes for lime in stats_types]
    header = ','.join(row_names+fp_names+m_names+tp_names+value_names+stats_names)
    with open(outfile,'w') as g:
        g.write(header+'\n')
        for record in sorted(motif_report.records.values(), key=lambda x: (x.transcript, x.query, x.start)):
            info,bases = record.generate_name().split('_'),list(record.generate_seq())
            A,B,D = record.A_react(),record.B_react(),record.D_react()
            metrics = ['change','increases','decreases','abs_change']
            stats = ['_'.join([prefix,stat]) for prefix in prefixes for stat in metrics]
            values = [getattr(record,stat) for stat in stats]
            numbs = [str(x) for x in A+B+D+values]
            g.write(','.join(info+bases+numbs)+'\n')
 
#Main Function
def main():
    parser = argparse.ArgumentParser(description='Searches for reactivity differences around given motifs')
    parser.add_argument('control',type=str,help='control <.react> file')
    parser.add_argument('experimental',type=str,help='experimental <.react> file')
    parser.add_argument('fasta',type=str,help='<.fasta> to pull sequences from')
    parser.add_argument('motif',type=str,help='Input file or motif')
    parser.add_argument('-fp',default=5,type=int, help='[default = 5] Bases to include 5\' of the motif')
    parser.add_argument('-tp',default=5,type=int, help='[default = 5] Bases to include 3\' of the motif')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-outdir',type=str,default = 'motif_out',help='[default = motif_out] Out Directory')
    parser.add_argument('-fastaout',action='store_true',default=False,help='Write windows to <.fasta> format as well')
    parser.add_argument('-reactout',action='store_true',default=False,help='Write accompanying <.react> files as well')
    args = parser.parse_args()
    
    #Check if path exists
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    #Read in reactivities and fasta
    control,experimental = map(structure_io.read_react,[args.control,args.experimental])
    sequences = structure_io.read_fasta(args.fasta)

    #Apply filter if input
    if args.restrict:
        covered = structure_io.read_restrict(args.restrict)
        sequences = {name:seq for name,seq in sequences.items() if name in covered}
    
    #Read in motifs
    motifs = read_motif(args.motif)
    
    #Out Nomenclature
    name_block_1 = [zz.replace('.react','') for zz in [args.control,args.experimental]]
    name_block_2 = [str(qq)+q for qq,q in zip([args.fp,args.tp],['fp','tp'])]
    
    #Iterate through motif(s)
    for motif in motifs:
        
        #Create a full MotifReport
        report = MotifReport(motif,sequences,control,experimental,args.fp,args.tp)
        
        #Generate name for each outfile
        out_name = '_'.join(name_block_1+[report.motif]+name_block_2)+'.csv'
        
        #Write out motif <.csv>
        write_motif_csv(report,os.path.join(args.outdir,out_name))

        #Write out motif <.fasta>
        if args.fastaout:
            out_fasta_name = out_name.replace('.csv','.fasta')
            fasta_dict = {r.generate_name():r.generate_seq() for r in report.records.values()}
            structure_io.write_fasta(fasta_dict,os.path.join(args.outdir,out_fasta_name))

        #Write out motif <.react>
        if args.reactout:
            control_out = {c.generate_name():c.A_react() for c in report.records.values()}
            exp_out = {e.generate_name():e.B_react() for e in report.records.values()}
            control_new = '_'.join([args.control.replace('.react',''),motif,str(args.fp)+'fp',str(args.tp)+'tp'])+'.react'
            exp_new = '_'.join([args.experimental.replace('.react',''),motif,str(args.fp)+'fp',str(args.tp)+'tp'])+'.react'
            structure_io.write_react(control_out,os.path.join(args.outdir,control_new))
            structure_io.write_react(exp_out,os.path.join(args.outdir,exp_new))

if __name__ == '__main__': 
    main()
