#!/usr/bin/env python2

#Imports
import argparse
from sf2libs import structure_io
from sf2libs.react_utils import generate_motif_coords_basic,remove_overlaps,subtract_reacts,delta_metrics
from sf2libs.basic_utils import flatten_list

#Classes
class MetaEntry(object):
    '''Just holds all the components of an entry in a clean way'''
    def __init__(self,transcript,coords,*categories):
        self.transcript = transcript
        self.pycoords = coords
        self.start = coords[0]+1
        self.end = coords[1]+1
        for sub_dict in categories:
            for k in sub_dict:
                setattr(self, k, sub_dict[k])

    def generate_key(self):
        return (self.transcript,self.start,self.end)

    def generate_name(self):
        return '_'.join([self.transcript,str(self.start),str(self.end)])

    def generate_seq(self):
        temp = self.fp_seq+self.m_seq+self.tp_seq
        return temp.replace('-','')

    def A_react(self):
        return filter(lambda R: R !='-',self.fp_A+self.A+self.tp_A)
 
    def B_react(self):
        return filter(lambda R: R !='-',self.fp_A+self.A+self.tp_A)

class MetaMotifReport(object):
    '''Generates and Contains a full MetaMotif Report'''
    def __init__(self,motifs,n_motifs,m_window,FP,TP,overlap,fasta_dict,react_A,react_B,fill='-'):
        self.motifs = motifs
        self.n_motifs = n_motifs
        self.m_window = m_window
        self.overlap = overlap
        self.FP = FP
        self.TP = TP
        self.records = {}
        #Populate Records
        for transcript,sequence in fasta_dict.items():
            try:
                A,B = react_A[transcript],react_B[transcript]
                delta = subtract_reacts(B,A)
                matches = scrampler(sequence,motifs,n_motifs,m_window,overlap)
                vectors = [A,B,delta,sequence]
                
                for m in matches:
                    meta_range = (m[0][0],m[-1][-1])
                    subseq = sequence[meta_range[0]:meta_range[1]]
                    gap_coords = [(m[i][1],m[i+1][0]) for i in range(0,len(m)-1)]
                    metas = [sequence[AA:Z] for AA,Z in m]
                    gaps = [sequence[a:z].lower() for a,z in gap_coords]
                    synth = ''.join([metas[i]+gaps[i] if i<len(gaps) else metas[i] for i in range(0,len(metas))])
                    mA,mB = A[meta_range[0]:meta_range[1]],B[meta_range[0]:meta_range[1]]
    
                    #Generate meta-motif sub-metrics
                    meta_values = [delta[BB:Y] for BB,Y in m]
                    gap_values = [delta[b:y] for b, y in gap_coords]
                    meta_metrics = delta_metrics(flatten_list(meta_values),'meta')
                    gap_metrics = delta_metrics(flatten_list(gap_values),'gap')
                    izaro_items = {'superseq':synth,'meta_coords':metas,'meta_chunks':metas,'m_seq':subseq,'A':mA,'B':mB}

                    #Generate 5'/3' stuff
                    fpA,fpB,fpD,fpS = [fp[max(meta_range[0]-FP,0):meta_range[0]] for fp in vectors]
                    tpA,tpB,tpD,tpS = [tp[meta_range[1]:min(len(tp)-1,meta_range[1]+TP)] for tp in vectors]

                    if len(fpS) < FP:
                        fpfill = (FP-len(fpS))*fill
                        fpS = fpfill+fpS
                        fpA,fpB,fpD = [list(fpfill)+orange for orange in [fpA,fpB,fpD]]

                    if len(tpS) < TP:
                        tpfill = (TP-len(tpS))*fill
                        tpS = tpS+tpfill
                        tpA,tpB,tpD = [lemon + list(tpfill) for lemon in [tpA,tpB,tpD]]

                    fp_items = {'fp_A':fpA,'fp_B':fpB,'fp_delta':fpD,'fp_seq':fpS}
                    tp_items = {'tp_A':tpA,'tp_B':tpB,'tp_delta':tpD,'tp_seq':tpS}
                    fp_items.update(delta_metrics(fpD,'fp'))
                    tp_items.update(delta_metrics(tpD,'tp'))

                    #Create Record, save
                    entry = MetaEntry(transcript,meta_range,fp_items,tp_items,izaro_items,meta_metrics,gap_metrics)
                    self.records[entry.generate_key()] = entry
            except KeyError:
                continue

#Functions
def scrampler(seq,motifs,nmotif=3,dist=50,allow_meta_overlaps=False):
    '''Ah. yes. The scrampler.'''
    raw_matches = sorted([z for m in motifs for z in generate_motif_coords_basic(seq,m)])
    matches = remove_overlaps(raw_matches)
    #
    collapsed = {}
    for sub_motif in matches:
        within_window = [x for x in matches if sub_motif[0] <x[0] < sub_motif[0]+dist]
        if len(within_window) >= nmotif:
            collapsed[(within_window[0][0],within_window[-1][-1])] = within_window
    #
    if allow_meta_overlaps:
        return sorted(collapsed.values())
    else:
        keep_keys = remove_overlaps(sorted(collapsed.keys()))
        hyper_collapsed = {k:v for k, v in collapsed.items() if k in keep_keys}
        return sorted(hyper_collapsed.values())

def write_meta_csv(motif_report,outfile='file.csv'):
    '''Writes out the default output'''
    row_names = ['transcript','metamotif','start','end','sequence']
    prefixes,stats_types = ['meta','gap','fp','tp'],['change','increases','decreases','delta']
    metrics = ['change','increases','decreases','abs_change']
    stats_names = [lemon+'_'+lime for lemon in prefixes for lime in stats_types]
    header = ','.join(row_names+stats_names)
    with open(outfile,'w') as g:
        g.write(header+'\n')
        for record in sorted(motif_report.records.values(), key=lambda x: (x.transcript, x.start)):
            base = [record.transcript,'|'.join(record.meta_chunks),str(record.start),str(record.end),record.superseq]
            stats = ['_'.join([prefix,stat]) for prefix in prefixes for stat in metrics]
            values = [getattr(record,stat) for stat in stats]
            numbs = [str(x) for x in values]
            g.write(','.join(base+numbs)+'\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Searches for reactivity differences around given multi-motifs')
    parser.add_argument('control',type=str,help='control <.react> file')
    parser.add_argument('experimental',type=str,help='experimental <.react> file')
    parser.add_argument('fasta',type=str,help='<.fasta> to pull sequences from')
    parser.add_argument('motifs',nargs='+',help='multi-motif components')
    parser.add_argument('-unique',action='store_false',default=True,help='[default = True] Remove overlapping multi-motifs')
    parser.add_argument('-mn',default=3,type=int,help='[default = 3] Number of multi-motifs required')
    parser.add_argument('-mw',default=50,type=int,help='[default = 50] Query window size for multi-motifs')
    parser.add_argument('-FP',default=30,type=int,help='[default = 30] Bases 5\' of multi-motifs')
    parser.add_argument('-TP',default=30,type=int,help='[default = 30] Bases 3\' of multi-motifs')
    parser.add_argument('-restrict',default=None,help='<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-fastaout',action='store_true',default=False,help='Write windows to <.fasta> format as well')
    parser.add_argument('-reactout',action='store_true',default=False,help='Write accompanying <.react> files as well')
    args = parser.parse_args()
    
    #Read in reactivities and fasta
    control,experimental = map(structure_io.read_react,[args.control,args.experimental])
    sequences = structure_io.read_fasta(args.fasta)

    #Apply filter if input
    if args.restrict:
        covered = structure_io.read_restrict(args.restrict)
        sequences = {name:seq for name,seq in sequences.items() if name in covered}
    
    #Out Nomenclature
    name_block_1 = [zz.replace('.react','') for zz in [args.control,args.experimental]]
    name_block_2 = [str(qq)+q for qq,q in zip([args.mn,args.mw,args.FP,args.TP],['mn','mw','FP','TP'])]
    
    #Generate base out name
    detail = 'shared' if args.unique else 'unique'
    out_name = '_'.join(name_block_1+name_block_2+[detail]+args.motifs)+'.csv'
    
    #Generate a full MetaReport
    report = MetaMotifReport(args.motifs,args.mn,args.mw,args.FP,args.TP,args.unique,sequences,control,experimental)

    #Write out motif <.csv>
    write_meta_csv(report,out_name)
    
    #Write out motif <.fasta>
    if args.fastaout:
        out_fasta_name = out_name.replace('.csv','.fasta')
        fasta_dict = {r.generate_name():r.generate_seq() for r in report.records.values()}
        structure_io.write_fasta(fasta_dict,out_fasta_name)

    if args.reactout:
        control_out = {c.generate_name():c.A_react() for c in report.records.values()}
        exp_out = {e.generate_name():e.B_react() for e in report.records.values()}
        #Writeout Control
        control_namae = [args.control.replace('.react','')]+name_block_2+[detail]+args.motifs
        control_new = '_'.join(control_namae)+'.react'
        structure_io.write_react(control_out,control_new)
        #Writeout Experimental
        exp_namae = [args.experimental.replace('.react','')]+name_block_2+[detail]+args.motifs
        exp_new = '_'.join(exp_namae)+'.react'
        structure_io.write_react(exp_out,exp_new)

if __name__ == '__main__':
    main()
