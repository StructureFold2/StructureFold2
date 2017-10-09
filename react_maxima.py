#!/usr/bin/python

'''
Takes two reacts - reports how many maxima are shared between them, i.e. there should be maxima changing if there
is a change in RNA structure between two reacts/conditions. If there are fewer total maxima in either <.react> than are
used to run, the result will be 'NA': impossible to share or not share what plain isn't there.
'''

#Imports
from itertools import islice
import argparse


#Functions
def read_in_target_transcripts(txt_file):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(txt_file,'r') as f:
        for line in f:
            good[line.strip()] = 'NULL'
    return good

def read_derived_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            reactivities = [float(x) if x != 'NA' else 'NA' for x in reactivities.split()]
            information[transcript] = reactivities
    return information

def maxima_indexer(alyst):
    '''Returns a sorted index of all the values in a list'''
    return sorted([z for z in zip(alyst, range(0,len(alyst))) if z[0] !='NA' and z[0] > 0.0 ],reverse=True)

def maxima_compare(alyst_1,alyst_2,overlap_keys,number_maxima=20,wiggle=3):
    '''does the thing'''
    data = {}
    
    if overlap_keys != None:
        for k in overlap_keys:
            try:
                a = set([x[1] for x in maxima_indexer(alyst_1[k])[0:number_maxima]])
                b = set([x[1] for x in maxima_indexer(alyst_2[k])[0:number_maxima]])
                if len(a) < number_maxima or len(b) < number_maxima:
                    data[k] = 'NA'
                    continue
                first_intersect =  a.intersection(b)
                #We want a little wiggle room
                wiggles = []
                for number in a.difference(b):
                    sub_wiggles = range(number-wiggle,number+1+wiggle)
                    wiggles = sub_wiggles+wiggles
                clean_wiggles = set(wiggles).difference(a)
                second_intersect = clean_wiggles.intersection(b)
                total_intersect = first_intersect.union(second_intersect)
                data[k] = list(total_intersect)
            except KeyError:
                print 'Entry for {} not found in one or more <.react> files, ignoring...'.format(k)
                continue
    else:
        #Get shared keys, run regardless
        q_keys = set(alyst_1.keys()).intersection(set(alyst_2.keys()))
        for k in q_keys:
            a = set([x[1] for x in maxima_indexer(alyst_1[k])[0:number_maxima]])
            b = set([x[1] for x in maxima_indexer(alyst_2[k])[0:number_maxima]])
            if len(a) < number_maxima or len(b) < number_maxima:
                data[k] = 'NA'
                continue
            first_intersect =  a.intersection(b)
            #We want a little wiggle room
            wiggles = []
            for number in a.difference(b):
                sub_wiggles = range(number-wiggle,number+1+wiggle)
                wiggles = sub_wiggles+wiggles
            clean_wiggles = set(wiggles).difference(a)
            second_intersect = clean_wiggles.intersection(b)
            total_intersect = first_intersect.union(second_intersect)
            data[k] = list(total_intersect)
    return data

def dump_csv(adict,outfile='shared_maxima.csv'):
    '''Writes a <.csv>'''
    with open(outfile,'w') as g:
        g.write(','.join(['transcript','shared_maxima'])+'\n')
        for k, v in adict.items():
            info_bit = str(len(v)) if v != 'NA' else 'NA'
            g.write(','.join([k,info_bit])+'\n')


def main():
    parser = argparse.ArgumentParser(description='Compares the reactivity maxima of every transcript between two <.react> files')
    parser.add_argument("react1",type=str,help="<.react> File for condition 1")
    parser.add_argument("react2",type=str,help="<.react> File for condition 2")
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-maxima',type=int,default=20, help='<int> [default = 20] Number of maxima to pick from both transcripts')
    parser.add_argument('-wiggle',type=int,default=3, help='<int> [default = 3] Number of bases maxima can be off by')
    args = parser.parse_args()
    
    #Read in a restrictive list, if provided
    restrict_dict = None if args.restrict == None else read_in_target_transcripts(args.restrict)
    
    #Read in <.rtsc>
    r_dict_1, r_dict_2 = read_derived_reactivities(args.react1),read_derived_reactivities(args.react2)
    
    #Cross Compare
    fused = maxima_compare(r_dict_1,r_dict_2,restrict_dict,args.maxima,args.wiggle)
    
    #Write out
    name = '_'.join(sorted([x.replace('.rtsc','') for x in [args.react1,args.react2]])+['maxima',str(args.maxima),str(args.wiggle)])+'.csv'
    dump_csv(fused,name)


if __name__ == '__main__':
    main()
    
    

    



