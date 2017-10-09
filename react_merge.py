#!/usr/bin/env python

'''
More Directions go here at some point
'''

#Imports
import argparse
from itertools import islice 

def read_in_react(react_fyle):
    '''Reads a <.react> file into a dictionary'''
    information = {}
    with open(react_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,stops = [n.strip() for n in next_n_lines]
            information[transcript] = [x for x in stops.split('\t')]
    print '\033[92mInput file\033[0m {} \033[92m contained\033[0m {}\033[92m entries\033[0m'.format(react_fyle,str(len(information)))
    return information

def reactivity_combine_max(react_1,react_2):
    '''Combines'''
    out_info = {}
    for transcript, values in react_1.items():
        if transcript in react_2:
            x_combine,x_new = zip(values,react_2[transcript]),[]
            for item in x_combine:
                sub_values = [float(x) for x in item if x != 'NA']
                if sub_values == []:
                    x_new.append('NA')
                else:
                    x_new.append(str(max(sub_values)))
            out_info[transcript] = x_new
        else:
            out_info[transcript] = values
    return out_info

def reactivity_combine_gaps(react_1,react_2):
    '''Combines, only replacing NAs in react_1 if there's a value in react_2'''
    out_info = {}
    for transcript, values in react_1.items():
        if transcript in react_2:
            x_combine,x_new = zip(values,react_2[transcript]),[]
            for item in x_combine:
                old,new = item[0], item[1]
                if old == 'NA' and new != 'NA':
                    x_new.append(new)
                else:
                    x_new.append(old)
            out_info[transcript] = x_new
        else:
            out_info[transcript] = values
    return out_info

def write_react(react_info,outfyle='out.react'):
    '''Writes the <.react> back out'''
    with open(outfyle,'w') as g:
        for transcript, data in react_info.items():
            g.write(transcript+'\n')
            g.write('\t'.join(data)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Combines <.react> files, takes the maximum value.')
    parser.add_argument('alpha',type=str,help='File 1, file you wish to supplement')
    parser.add_argument('beta',type=str,help='File 2, file to be added to File 1')
    parser.add_argument('-gaps',action='store_true',default=False,help='Only fill in NAs of file 1 with file 2')
    parser.add_argument('-name',default=None, help='Specify output file name')
    args = parser.parse_args()
    
    #Generate out name
    default_name = '~'.join([args.alpha.replace('.react',''),args.beta.replace('.react','')])+'.react'
    out_name = default_name if args.name == None else args.name
    
    #Read in data to be combined
    alpha_dict,beta_dict = read_in_react(args.alpha),read_in_react(args.beta)
    
    #Combine Data
    new_reacts = reactivity_combine_max(alpha_dict,beta_dict) if args.gaps == False else reactivity_combine_gaps(alpha_dict,beta_dict)
    
    #Write Data
    write_react(new_reacts,out_name)


if __name__ == '__main__':
    main()
