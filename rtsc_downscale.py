#!/usr/bin/env python2

'''
Downscales <.rtsc> files, keeping a given percentage of the original RT hits.
This downscaling cab be done statically, simply multiplying every RT count by a specified percentage
This downscaling can be done pseudo-randomly, removing hits randomly down to a specified percentage.
This downscaling can be done pseudo-randomly, continually checking a random bases and subtracting a 
stop until enough have been removed.
'''

#Imports
from sf2libs.structure_io import read_rtsc, write_rtsc, read_restrict
import random
import glob
import argparse

#Functions
def keep_static_percentage(alist,percentage=.50):
    '''Evenly retains a given percentage of RT hits. Can generate fractional hits.'''
    return [float(i)*percentage for i in alist]

def keep_random_downsample(num_list,percentage=.50):
    '''Randomly retains a given percentage of RT hits'''
    total,length,internal = sum(num_list),len(num_list),num_list[:]
    keeps = int(total*percentage)
    removes = total-keeps
    while removes > 0:
        target = random.randrange(length)
        if internal[target] > 0:
            internal[target]-=1
            removes+=-1
    return internal

def random_per_stop(alist,chance=.50):
    '''Roll a Die for each RT stop at each positon'''
    newlist,new_value = [],int(chance*1000)
    for value in alist:
        dice_rolls = [random.randint(1,1000) for stop in range(0,value)]
        keeps = filter(lambda x:x <= new_value, dice_rolls)
        newlist.append(len(keeps))
    return newlist

def main():
    parser = argparse.ArgumentParser(description='Downscales <.rtsc> files.')
    parser.add_argument('-f',default=None, nargs='+', help = 'Specific <.rtsc> to operate on')
    parser.add_argument('mode',type=str.upper, choices = ['FRACTIONAL','RANDOMREAD','RANDOMPOSITION'])
    parser.add_argument('-ratio',type=float,default=.50, help='[default = 0.50] Fraction of RT stops to retain')
    parser.add_argument('-restrict',default = None, help = 'Limit downscaling to these specific transcripts <.txt> ')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    args = parser.parse_args()

    #Files to operate on, dictionary of functions
    fyle_lyst = sorted(glob.glob('*.rtsc')) if args.f == None else sorted(args.f)
    downscale_methods = {'FRACTIONAL':keep_static_percentage,'RANDOMPOSITION':keep_random_downsample,'RANDOMREAD':random_per_stop}

    #Iterate through file(s) 
    for fyle in fyle_lyst:
        
        #Read in the <.rtsc>, generate new name
        data = read_rtsc(fyle)
        new_fyle = fyle.replace('.rtsc','_'+args.mode+'_'+str(args.ratio).replace('.','')+'.rtsc')
        
        #You could be doing random for a long time if you do not filter by coverage.
        if args.restrict != None:
            covered = read_restrict(args.restrict)
            data = {name:stops for name,stops in data.items() if name in covered} 

        new_data = {k:downscale_methods[args.mode](v,args.ratio) for k, v in data.items()}
        write_rtsc(new_data,new_fyle,args.sort)

#Main
if __name__ == '__main__':
    main()
