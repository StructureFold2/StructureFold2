#!/usr/bin/python

'''
Converts a <.rtsc> to a <.csv> of relative transcript abundances (RTSC hits per kilobase per million reads). 
These abundances are approximate and traditional RNA-seq is more precise at this task than Structure-Seq libraries, 
but for a quick analysis or inlieu of a parallel RNA-seq experiment, this tool will yeild a good approximation.
Current math:
Number of RSTC hits*1000*1000000/total reads in file*length of transcript
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
            information[transcript] = (sum(int(x) if x !='NA' else 0.0 for x in stops.split('\t')),len(stops))
    return information

def dump_dictionary(adict,outfyle,total_n):
    '''Writes out the data'''
    with open(outfyle,'w') as g:
        g.write(','.join(['transcript',outfyle.replace('_abundances.csv','')+'_'+'RTPKM'])+'\n')
        for k, v in adict.items():
            number = (float(v[0])*1000*1000000)/(v[1]*total_n)
            g.write(','.join([k,str(number)])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Determines approximate transcript abundance based on <.rtsc> files.')
    parser.add_argument('-rtsc',default = None, help = 'Operate on these files, rather than the directory')
    args = parser.parse_args()
    #Files to operate on
    fyle_lyst = sorted(glob.glob('*.rtsc')) if args.rtsc == None else sorted(args.rtsc)
    #Loop
    for fyle in fyle_lyst:
        data = read_in_total_stops(fyle)
        total_number = sum(v[0] for v in data.values())
        new_fyle = fyle.replace('.rtsc','')+'_abundances.csv'
        dump_dictionary(data,new_fyle,total_number)

if __name__ == '__main__':
    main()
