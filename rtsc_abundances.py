#!/usr/bin/env python

'''
Converts <.rtsc> to <.csv> of relative transcript abundances. 
Traditional RNA-seq may be more precise at this task than Structure-Seq libraries.
'''

#Imports
from itertools import islice
import glob
import argparse

def read_in_total_stops(afile):
    '''Reads in <.rtsc>, returns the total stops and length for each transcript'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,derp = [n.strip() for n in next_n_lines]
            expanded_stops = [float(x) for x in stops.split('\t')]
            information[transcript] = (sum(expanded_stops),len(expanded_stops))
    return information

def values_to_RPKM(data):
    '''Calculates RPKM or Reads Per Kilobase per Million reads'''
    normalize,RPKM_values = sum(v[0] for v in data.values()),{}
    for k, v in data.items():
        RPKM_value = (float(v[0])*1000*1000000)/(v[1]*normalize)
        RPKM_values[k] = RPKM_value
    return RPKM_values

def values_to_TPM(data):
    '''Calculates TPM or Transcripts Per Kilobase Million reads'''
    TPM_values = {}
    reads_per_kb = dict([(k,v[0]/(float(v[1])/1000)) for k, v in data.items()])
    normalize = sum(reads_per_kb.values())/1000000
    for k, v in reads_per_kb.items():
        TPM_values[k] = v/normalize
    return TPM_values

def write_data(adict,outfyle,data_unit):
    '''Writes the data to a <.csv>'''
    header = ','.join(['transcript',outfyle.replace('.csv','')])
    with open(outfyle,'w') as g:
        g.write(header+'\n')
        for transcript,abundance_stat in adict.items():
            g.write(','.join([transcript,str(abundance_stat)])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Determines approximate transcript abundance based on <.rtsc> files.')
    parser.add_argument('mode',type=str.upper, choices = ['RPKM','TPM'])
    parser.add_argument('-f',default = None, nargs='+', help = 'Specifc <.rtsc> to operate on')
    args = parser.parse_args()

    #Files to operate on, dictionary of functions
    fyle_lyst = sorted(glob.glob('*.rtsc')) if args.f == None else sorted(args.f)
    abundance_methods = {'RPKM':values_to_RPKM,'TPM':values_to_TPM}

    #Iterate through file(s) 
    for fyle in fyle_lyst:
        data,new_fyle = read_in_total_stops(fyle),fyle.replace('.rtsc','_'+args.mode+'.csv')
        processed_data = abundance_methods[args.mode](data)
        write_data(processed_data,new_fyle,args.mode)


if __name__ == '__main__':
    main()






 
