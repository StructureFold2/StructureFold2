#!/usr/bin/env python

'''
This script is to be run on the output CT directory (or any other directory full of <.ct> files) of a batch fold
where the -MFE setting was not called; the MFE structure will be the first reported structure, and thus be the first
entry in the <.ct> file. This script will grab the free energy of all MFE structures and write them to a <.csv>.
'''

#Imports
import os
import subprocess
import re
import argparse

#Functions
def glean_ct_tops(directory):
    '''Grabs the top line of every file in the <.ct> directory, returns dictionary'''
    home,info = os.getcwd(),{}
    os.chdir(directory)
    procced_info = subprocess.Popen('head -1 *.ct', shell=True, stdin=subprocess.PIPE, 
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True).stdout.read()
    read_data = [line for line in procced_info.split('\n') if line != '']
    for i in range(0,len(read_data),2):
        pair = read_data[i:i+2]
        key = '_'.join(pair[0].strip().split()[1].split('_')[:-2])
        try:
            value = re.findall("-\d+\.\d+", pair[1])[0]
        except IndexError:
            value = 0
        info[key] = value
    os.chdir(home)
    return info

def write_csv(adict,outfyle):
    '''Writes out a simple <.csv> with the MFE data'''
    with open(outfyle, 'w') as out:
        out.write(','.join(['transcript','Free_Energy'])+'\n')
        for k, v in adict.items():
            out.write(','.join([k,str(v)])+'\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Takes a CT directory of a batch fold, writes a <.csv> of transcript MFEs.')
    parser.add_argument('directory',default = None, help = 'Operate on this directory')
    parser.add_argument('-name',default = 'MFE.csv', help = '[default = MFE.csv] Output file name')
    args = parser.parse_args()
    #
    MFEdict = glean_ct_tops(args.directory)
    write_csv(MFEdict,args.name)


if __name__ == '__main__': 
    main()
