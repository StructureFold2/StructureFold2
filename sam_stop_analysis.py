#!/usr/bin/env python
#David Tack

#Imports
import collections
import argparse
import itertools
from sf2libs.structure_io import read_fasta,check_extension

class RT_Report(object):
    '''An RT report for a transcript'''
    def __init__(self,name,sequence,stops=None,prior=None):
        self.name = name
        self.sequence = sequence
        self.rt_stops = collections.Counter() if not stops else stops
        self.rt_prior = collections.Counter() if not prior else prior

    def add_transcript(self,position,read):
        '''Adds the read/RT stop to the prior/stops'''
        new_pattern = collections.Counter([read[n-1:n+1] for n in xrange(1,len(read))])
        self.rt_prior = self.rt_prior + new_pattern
        if position != 1:
            stop = self.sequence[position-2:position]
            self.rt_stops[stop]+=1

    def __str__(self):
        truncated_seq = self.sequence[0:5]+'...'+self.sequence[-5:]
        return '\n'.join([self.name,truncated_seq])
    
    def generate_rows(self,patterns):
        '''Summarize Stats'''
        p_sum = sum(self.rt_prior.values())
        s_sum = sum(self.rt_stops.values())
        p_vals = [self.rt_prior[a] for a in patterns]
        s_vals = [self.rt_stops[b] for b in patterns]
        p_fracts = [float(self.rt_prior[c])/p_sum if p_sum else 'NA' for c in patterns]
        s_fracts = [float(self.rt_stops[d])/s_sum if s_sum else 'NA' for d in patterns]
        enrichment = [s - p if all(isinstance(q,float) for q in [s,p]) else 'NA' for s, p in zip(s_fracts,p_fracts)]
        enrichment2 = [round(z*100,5) for z in enrichment]
        p_fracts2 = [round(e*100,5) for e in p_fracts]
        s_fracts2 = [round(f*100,5) for f in s_fracts]
        stats = zip(patterns,p_vals,p_fracts2,s_vals,s_fracts2,enrichment2)
        #return {stat[0]:stat[1:] for stat in stats}
        return stats

class Mapped_Read(object):
    '''Just holds a mapped read'''
    def __init__(self,a_sam_line):
        items = a_sam_line.strip().split('\t')
        self.qname = items[0]
        self.flag = int(items[1])
        self.rname = items[2]
        self.pos = int(items[3])
        self.mapq = int(items[4])
        self.cigar = items[5]
        self.rnext = items[6]
        self.pnext = items[7]
        self.tlen = items[8]
        self.seq = items[9]
        self.qual = items[10]
        if len(items) > 11:
            self.attributes = items[11:] 

def populate_reports(empty_reports,files):
    '''Populates Reports'''
    for fyle in files:
        with open(fyle,'r') as temp:
            for line in temp:
                data = Mapped_Read(line)
                empty_reports[data.rname].add_transcript(position=data.pos,read=data.seq)

def generate_summary(rt_reports):
    '''Generates an amalgamated entry'''
    all_stops, all_prior = collections.Counter(),collections.Counter()
    for entry in rt_reports.values():
        all_stops = all_stops + entry.rt_stops
        all_prior = all_prior + entry.rt_prior
    return RT_Report(name='Summary',sequence='ABCDEF',stops=all_stops,prior=all_prior)

def write_report(data,patterns,threshold,name,des):
    '''Herpderp'''
    with open(name,'w') as g:
        col_labels = ['prior','priorperc','stops','stopsperc','enrichment']
        col_names = ['_'.join([des,label]) for label in col_labels] if des else col_label
        g.write(','.join(['transcript','pattern']+col_names)+'\n')
        for entry in data.values():
            if sum(entry.rt_stops.values()) > threshold:
                rows = entry.generate_rows(patterns)
                for row in rows:
                    out_row = ','.join([entry.name]+[str(i) for i in row])
                    g.write(out_row+'\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Transcript-Specific RT stop patterns')
    parser.add_argument('index',type=str,help='Reference Fasta File')
    parser.add_argument('-f',default = None, nargs='+', help = 'Specifc <.sam> to operate on')
    parser.add_argument('-desc',type=str,default = '', help='Column Header Prefix')
    parser.add_argument('-total',action='store_true', help='Calculate a summary entry')
    parser.add_argument('-thres',type=int, default = 100, help = 'Required number of RT stops to report')
    parser.add_argument('-name',type=str,default = 'outfile.csv', help='[default = outfile.csv] Name of the output')
    args = parser.parse_args()
    
    cannonical = [''.join(q) for q in itertools.product('ACGT',repeat=2)]
    
    #Workflow
    out_name = check_extension(args.name,'.csv')
    
    #Read in reference transcriptome
    fasta_data = read_fasta(args.index)
    print('Found {} sequence records in {}'.format(len(fasta_data),args.index))
    
    #Generate all potential RT reports
    rt_data = {k:RT_Report(name=k,sequence=v) for k, v in fasta_data.items()}
    print('Populated {} transcript reports').format(len(rt_data))
    
    #Populate RT reports
    populate_reports(empty_reports=rt_data,files=args.f)
    
    #Generate summary report
    if args.total:
        summary_data = {'Summary':generate_summary(rt_data)}
        summary_name = out_name.replace('.csv','_summary.csv')
        write_report(summary_data,cannonical,args.thres,summary_name,args.desc)

    #Write Out
    transcripts_name = out_name.replace('.csv','_transcripts.csv')
    write_report(rt_data,cannonical,args.thres,transcripts_name,args.desc)

if __name__ == '__main__': 
    main()