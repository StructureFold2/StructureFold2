#!/usr/bin/env python

'''
A simple class for working with, and extracting data from, Connectivity Tables <.ct> files
Contains the CT_Report class, as well as collect_connectivity_tables(directory) to collect your ct data
'''

#Imports
import os
import glob

class CT_Report(object):
    '''Takes in CT files, make them useable for a bit.'''
    def __init__(self,a_ct_file):
        with open(a_ct_file,'r') as f:
            header_bits = f.next().strip().split()
            self.length = int(header_bits[0])
            self.name = header_bits[1]
            lines = [line.strip().split() for line in f]
            n,base,bl,br,paired,n2 = [[q[i] for q in lines] for i in range(0,6)]
            self.flat_index = n
            self.base_index = base
            self.left = bl
            self.right = br
            self.paired_index = paired
            self.para_index = n2
            self.perc_stranded = len([z for z in paired if z != '0'])/float(int(header_bits[0]))
            self.perc_unstranded = paired.count('0')/float(int(header_bits[0]))
    
    def primary_strand(self,coord_tuple):
        '''Takes Pythonic coordinates, returns the nucleotide sequence of the primary strand'''
        try:
            return ''.join(self.base_index[coord_tuple[0]:coord_tuple[1]])
        except IndexError:
            return None
    
    def complement_strand(self,coord_tuple):
        '''Takes Pythonic coordinates, returns the nucleotide sequence of any complement to the fold'''
        try:
            pairings = [int(z) for z in self.paired_index[coord_tuple[0]:coord_tuple[1]]]
            paired_bases = [self.base_index[g-1] if g > 0 else '-' for g in pairings]
            #The -1 is because the file is non Pythonic
            return ''.join(paired_bases)
        except IndexError:
            return None
    
    def complement_loops(self):
        '''Returns a list of any unpaired stretches of nucleotides from the complementary strand'''
        paired_bases = [self.base_index[g-1] if g > 0 else '-' for g in [int(z) for z in self.paired_index]]
        dummy_bases = ''.join([q if q == '-' else '.' for q in paired_bases]).split('.')
        cleaned_loops = [z for z in dummy_bases if z != '']
        return cleaned_loops

    def complement_stems(self):
        '''Returns a list of any paired stretches of nucleotides from the complementary strand'''
        paired_bases = [self.base_index[g-1] if g > 0 else '-' for g in [int(z) for z in self.paired_index]]
        cleaned_stems = [z for z in ''.join(paired_bases).split('-') if z != '']
        return cleaned_stems

    def primary_loops(self):
        '''Returns a list of any unpaired stretches of nucleotides from the primary strand'''
        paired_bases = [self.base_index[g-1] if g > 0 else '-' for g in [int(z) for z in self.paired_index]]
        unpaired_bases = ['-' if paired_bases[i] != '-' else self.base_index[i] for i in range(0,len(paired_bases))]
        cleaned_up = [z for z in ''.join(unpaired_bases).split('-') if z != '']
        return cleaned_up

    def primary_stems(self):
        '''Returns a list of any paired stretches of nucleotides from the primary strand'''
        paired_bases = [self.base_index[g-1] if g > 0 else '-' for g in [int(z) for z in self.paired_index]]
        base_paired_bases = [self.base_index[i] if paired_bases[i] != '-' else '-' for i in range(0,len(paired_bases))]
        cleaned_up = [z for z in ''.join(base_paired_bases).split('-') if z != '']
        return cleaned_up

def collect_connectivity_tables(directory):
    '''Snags a directory worth of CT files'''
    home,data = os.getcwd(),{}
    os.chdir(directory)
    for fyle in glob.glob('*.ct'):
        data[fyle] = CT_Report(fyle)
    os.chdir(home)
    return data

def w_index_to_w_range(value,w_len,w_step,t_len=1000):
    '''Converts a step index back to a coordinate range, i.e. to pull features from react_windows'''
    blank_windows = [(i,i+w_len) for i in xrange(0,t_len-(w_len-1), w_step)]
    return blank_windows[value]

def main():
    pass

if __name__ == '__main__':
    main()
