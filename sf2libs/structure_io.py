#All of the file io stuff for StructureFold2

import itertools
from Bio import SeqIO

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string 

def read_fasta(afasta):
    '''Fasta to Python dictionary'''
    fasta_sequences = SeqIO.parse(open(afasta),'fasta')
    return {record.id:str(record.seq) for record in fasta_sequences}

def read_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary, transcript_name:[list of stop numbers]'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(itertools.islice(f,3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information

def read_react(react_file):
    '''Reads a <.react> file into a dictionary, transcript_name:[list of reactivities] '''
    information = {}
    with open(react_file,'r') as f:
        while True:
            next_n_lines = list(itertools.islice(f,2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!= 'NA' else 'NA' for x in reactivities.split()]
    return information

def read_restrict(coverage_overlap):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_overlap,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def write_fasta(info,outfyle='out.fasta',sort_flag=False,line_width=80):
    '''Writes out a dictionary as a <.fasta>, line_width controls chars per line'''
    with open(outfyle,'w') as g:
        if sort_flag:
            for name,seq in sorted(info.items()):
                g.write('>'+ name+'\n')
                for i in xrange(0,len(seq),line_width):
                    g.write(seq[i:i+line_width]+'\n')
        else:
            for name,seq in info.items():
                g.write('>'+ name+'\n')
                for i in xrange(0,len(seq),line_width):
                    g.write(seq[i:i+line_width]+'\n')

def write_react(react_dictionary,outfile='data.react',sort_flag=False):
    '''Writes out a dictionary as a <.react> file'''
    with open(outfile,'w') as g:
        if sort_flag:
            for transcript, entry in sorted(react_dictionary.items()):
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in entry])+'\n')
        else:
            for transcript, entry in react_dictionary.items():
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in entry])+'\n')

def write_rtsc(rtsc_dictionary,outfile='data.rtsc',sort_flag=False):
    '''Writes out a dictionary as a <.rtsc> file'''
    with open(outfile,'w') as g:
        if sort_flag:
            for transcript, entry in sorted(rtsc_dictionary.items()):
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in entry])+'\n\n')
        else:
            for transcript, entry in rtsc_dictionary.items():
                g.write(transcript+'\n')
                g.write('\t'.join([str(number) for number in entry])+'\n\n')


