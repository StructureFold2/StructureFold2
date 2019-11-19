#!/usr/bin/env python2

'''
Reformats <.react> files such that correlation may be easily caluclated, 
either genome-wide or on a per transcript basis.Output is a <.csv>.
You may filter output by using a coverage overlap, 
where you only want to get the correlation between transcripts mutually above a certain threshold.
'''

#Imports
import argparse
import itertools

#Functions
def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = 'NULL'
    return good

def read_in_react(react_fyle):
    '''Reads a <.react> file into a dictionary'''
    information = {}
    with open(react_fyle,'r') as f:
        while True:
            next_n_lines = list(itertools.islice(f, 2))
            if not next_n_lines:
                break
            transcript,stops = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x != 'NA' else 'NA' for x in stops.split()]
    print 'File: {} had a total of {} entries'.format(react_fyle,len(information))
    return information

def read_in_reacts(reacts):
    '''Reads in any number of reacts, derps repeatedly'''
    react_data,new = map(read_in_react,reacts),{}
    #Generate common keys among all files, only work with these.
    common_keys = set.intersection(*map(set, react_data))
    print 'Found {} common keys (transcripts) between all files'.format(len(common_keys))
    #Fill nested dictionary.
    for name, subdict in zip(reacts,react_data):
        for transcript, react_pattern in subdict.items():
            new.setdefault(transcript, {})[name.replace('.react','')] = react_pattern
    return new

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def write_react_repeatability(react_data,reacts_in,out_fyle):
   '''Writes out the columns'''
   sub_keys = [q.replace('.react','') for q in reacts_in]
   header = ','.join(['transcript']+sub_keys)
   with open(out_fyle,'w') as g:
       g.write(header+'\n')
       for transcript, data in react_data.items():
           #Assumption is that, if there's a missing transcript in the sub_dict
           #We will fill it in with up to 100k NAs. Note the stdout should alert the user to a discrepancy
           order = [data[key] if key in data else ['NA']*99999 for key in sub_keys]
           for i, mini in enumerate(zip(*order),1):
               g.write(','.join([transcript,str(i)]+[str(q) for q in mini])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Reformats <.react> for easy correlation analysis')
    parser.add_argument('react',help='Input <.react> files',nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-restrict',default=None, help='Filter to these transcripts via coverage file')
    args = parser.parse_args()
    
    #Read in data, filter out transcripts if applicable
    data = read_in_reacts(args.react)
    
    #Restrictions, if they exist
    restrict_dict = None if args.restrict == None else read_in_target_transcripts(args.restrict)
    if restrict_dict:
        data = dict([(x,z) for x,z in data.items() if x in restrict_dict])
    
    #Outfile nomenclature
    default_name = '_'.join([x.replace('.react','') for x in args.react])+'_react_correlation.csv'
    out_name = default_name if args.name == None else check_extension(args.name,'.csv')
    
    #Write Out
    write_react_repeatability(data,args.react,out_name)

if __name__ == '__main__': 
    main()
 