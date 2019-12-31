#!/usr/bin/env python2

'''
Generates numerical summaries of batch MFE foldings. Depending on the mode it either...
F)Summarizes multiple directories of MFE <.ct> files into a single <.csv> sorted by transcript.
R)Summarizes a directory of MFE <.ct> files into a single <.csv> sorted by file name.
P)Runs 'scorer' between multiple directories of MFE <.ct> files, generating a pairwise PPV file.
'''

#Imports
import os
import argparse
import itertools
import subprocess
from sf2libs.connectivity_table import collect_connectivity_tables
from sf2libs.structure_io import check_extension

#Functions
def fuse_connectivity_dirs(list_of_directories,offset,empty):
    '''Applys the CT_Report table to directories, gathers information'''
    mass_data = {}
    for directory in list_of_directories:
        raw_dict = collect_connectivity_tables(directory,empty)
        #We want to collapse entries around transcript between directories
        #Some transcript nomenclatures have underscores in them, thus an offset can work around it
        melted_dict = {'_'.join(k.split('_')[:1+offset]):v for k, v in raw_dict.items()}
        mass_data[directory.strip(os.sep).split(os.sep)[-1]] = melted_dict
    return mass_data

def write_out_fused(f_data,out_fyle='derp.csv'):
    '''Writes out Fused data'''
    vertical_keys = sorted(list(set.union(*map(set, f_data.values()))))
    horizontal_keys = sorted(f_data.keys())
    cols_type = ['double_stranded','single_stranded','DeltaG']
    cols_labels = [item+'_'+affix for item in horizontal_keys for affix in cols_type]
    with open(out_fyle,'w') as g:
        header = ','.join(['transcript']+cols_labels)
        g.write(header+'\n')
        for v in vertical_keys:
            garble_line = [[f_data[h][v].perc_stranded,f_data[h][v].perc_unstranded,f_data[h][v].energy] if v in f_data[h] else ['NA']*3 for h in horizontal_keys]
            flat_line = [str(item) for sublist in garble_line for item in sublist]
            outline = ','.join([v]+flat_line)
            g.write(outline+'\n')

def write_out_raw(raw_data,out_fyle='derp.csv'):
    '''Writes out raw data'''
    header = ','.join(['File','double_stranded','single_stranded','DeltaG'])
    with open(out_fyle,'w') as g:
        g.write(header+'\n')
        for k, v in raw_data.items():
            aline = ','.join([str(x) for x in (k, v.perc_stranded,v.perc_unstranded,v.energy)])
            g.write(aline+'\n')

def program_check(prog):
    '''Check if the thing is in the path and a file'''
    return any(os.access(os.path.join(path,prog),os.X_OK) and os.path.isfile(os.path.join(path,prog)) for path in os.environ["PATH"].split(os.pathsep))

def read_ppv(ppv_file):
    '''Reads a PPV file'''
    fyle_name = ppv_file.split(os.sep)[-1].strip('.txt')
    accepted,predicted,transcript = fyle_name.split('~')
    with open(ppv_file,'r') as f:
        lines = f.readlines()
        ppv = lines[4].split('=')[-1].strip().strip('%')
    return transcript, {'~'.join([accepted,predicted]):ppv}

def wrangle_ppvs(directories,offset,temp_dyr='temp_ppv'):
    '''Builds a nested dictionary of all the possible PPVs between transcripts between directories'''

    #Create a directory, all the PPV outputs will go in it.
    os.mkdir(temp_dyr)

    #Generate all the PPV files
    dir_labels = [directory.strip(os.sep).split(os.sep)[-1] for directory in directories]
    dir_info = {label:dyr for label,dyr in zip(dir_labels,directories)}
    dir_combos = itertools.permutations(dir_labels,2)
    for combo in dir_combos:
        lemon = {'_'.join(m.split('_')[:1+offset]):m for m in os.listdir(dir_info[combo[0]])}
        lemon = { yellow:citrus for yellow, citrus in lemon.items() if citrus.endswith('.ct')}
        orang = {'_'.join(m.split('_')[:1+offset]):m for m in os.listdir(dir_info[combo[1]])}
        orang = { orange:sour for orange, sour in orang.items() if sour.endswith('.ct')}
        common_keys = set.intersection(*map(set, [lemon,orang]))
        for transcript in common_keys:
            lemon_path = os.path.join(dir_info[combo[0]],lemon[transcript])
            orang_path = os.path.join(dir_info[combo[1]],orang[transcript])
            citrus = os.path.join(temp_dyr, '~'.join(list(combo)+[transcript])+'.txt')
            command = ' '.join(['scorer',orang_path,lemon_path,citrus])
            subprocess.call(command,shell=True)

    #Collect the data
    collected_data = {}
    for fyle in os.listdir(temp_dyr):
        seq, data = read_ppv(os.path.join(temp_dyr,fyle))
        collected_data.setdefault(seq,{}).update(data)

    #Remove the directory, since we have read in all the data
    subprocess.call(' '.join(['rm','-r','temp_ppv']),shell=True)
    return collected_data

def write_out_ppv(ppv_data,fyle='out_ppv.csv'):
    '''Writes out a consolidated PPV <.csv>'''
    keys = sorted(list(set.union(*map(set, ppv_data.values()))))
    header= ','.join(['transcript']+[k+'_ppv' for k in keys])
    with open(fyle,'w') as g:
        g.write(header+'\n')
        for transcript, sub in sorted(ppv_data.items()):
           line=','.join([transcript]+[sub[key] if key in sub else 'NA' for key in keys])
           g.write(line+'\n')

def main():
    parser = argparse.ArgumentParser(description='Summarizes or compares MFE <.ct> files')
    parser.add_argument('-d',type=str,help='CT directory/directories', nargs='+')
    parser.add_argument('-mode',type=str.upper,default = None,choices = ['F','R','P'],help='Fused/Raw/PPV statistics')
    parser.add_argument('-name',type=str,default = None, help = 'Output file name')
    parser.add_argument('-na',type=str,default = 'NA', help = '[default = NA] Null deltaG value')
    parser.add_argument('-offset',type=int,default = 0, help = 'Number of Underscores in Transcript Names')
    args = parser.parse_args()
    
    if args.mode == 'R':
        directory = args.d[0]
        if len(args.d) > 1:
            print('Raw mode can only process a single directory, processing {}'.format(directory))
        data = collect_connectivity_tables(directory,args.na)
        melted_name = directory.strip(os.sep).split(os.sep)[-1]
        default_name = '_'.join([melted_name,'statisics'])+'.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.csv')
        write_out_raw(data,out_name)

    if args.mode == 'F':
        data = fuse_connectivity_dirs(args.d,args.offset,args.na)
        default_name = '_'.join(sorted(data.keys())+['statistics'])+'.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.csv')
        write_out_fused(data,out_name)

    if args.mode == 'P':
        if program_check('scorer'):
            data = wrangle_ppvs(args.d,args.offset)
            ppv_dirs = [dyr.strip(os.sep).split(os.sep)[-1] for dyr in args.d]
            default_name = '_'.join(sorted(ppv_dirs)+['PPV'])+'.csv'
            out_name = default_name if args.name == None else check_extension(args.name,'.csv')
            write_out_ppv(data,out_name)
        else:
            print 'Unable to locate an excutable version of scorer.'

    if args.mode == None:
        print 'No mode selected!'

if __name__ == '__main__':
    main()
