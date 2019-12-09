import os
import glob
import re

class CT_Report(object):
    '''Reads in a <.ct> file, holds the information'''
    def __init__(self,ct_file,DeltaG='NA'):
        with open(ct_file,'r') as f:
            header = f.next()
            header_bits = header.strip().split()
            self.length = int(header_bits[0])
            self.name = header_bits[-1]
            energy_info = re.findall(r"ENERGY = [-+]?\d*\.\d+",header)
            self.energy = float(energy_info[0].split()[2]) if energy_info else DeltaG
            lines = [line.strip().split() for line in f]
            n,base,bl,br,paired,n2 = [[q[i] for q in lines] for i in range(0,6)]
            self.flat_index = n
            self.base_index = base
            self.left = bl
            self.right = br
            self.paired_index = paired
            self.para_index = n2
            #Check that the length reported is correct, if so use len of an index list thing.
            self.perc_stranded = len([z for z in paired if z != '0'])/float(self.length)
            self.perc_unstranded = paired.count('0')/float(self.length)

def collect_connectivity_tables(directory,empty_dg):
    '''Snags a directory worth of CT files'''
    home,data = os.getcwd(),{}
    os.chdir(directory)
    for fyle in glob.glob('*.ct'):
        data[fyle] = CT_Report(fyle,empty_dg)
    os.chdir(home)
    return data