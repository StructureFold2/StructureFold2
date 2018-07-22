#!/usr/bin/env python
#StructureFold2 batch running script for Bowtie2
#This script may be easily modified to run your read mapper of choice.
#This is an accessory script and is not a core part of StructureFold2, but you may find it useful.

#Imports
import glob
import subprocess
import time
import argparse
import re

#Functions
def bowtie_batcher(infile,mapping_index,threads_number,out_suffix,phred_adjust,multimap,fail_reads,log_flag):
    '''Runs the bowtie2 command'''
    outfile = infile.split('.')[0]+'_'+out_suffix+'.sam'
    index_bin,threads_bin = ['-x',mapping_index],['-p',str(threads_number)]
    phred_bin = [] if phred_adjust == False else ['--phred64']
    map_bin = ['-a'] if multimap == True else []
    file_bin = ['-q',infile, '>',outfile]
    fail_bin = [] if fail_reads == True else ['--no-unal']
    command = ' '.join(['bowtie2']+map_bin+threads_bin+fail_bin+phred_bin+index_bin+file_bin)
    if log_flag == False:
        subprocess.call(command,shell=True)
    if log_flag == True:
            procced_info = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, 
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
            return procced_info.stdout.read(),outfile,command
        
def write_log_file(info_bucket,outfile):
    '''Takes the info collected from stout and makes a nice consistent <.csv>'''
    sub_buckets = [info_bucket[x:x+6] for x in range(0,len(info_bucket),6)]
    x_header = ','.join(['fastq','sam','in_reads','mapped_fract','unique_reads','unique_fract','multi_reads',
                         'multi_fract','unmapped_reads','unmapped_fract','time_start','time_end','bowtie2_command'])
    with open(outfile,'w') as g:
        g.write(x_header+'\n')
        for rec in sub_buckets:
            num = re.findall(r'[-+]?\d*\.\d+|\d+',rec[4])
            out_fields = ','.join([rec[1],rec[2],num[0],num[12],num[6],num[7],num[9],
                                   num[10],num[3],num[4],rec[0],rec[5],rec[3]])
            g.write(out_fields+'\n')

def main():
    parser = argparse.ArgumentParser(description='Runs Bowtie2 on every <.fastq> in the current directory using the recommended settings for a Structure-seq analysis.')
    parser.add_argument("index",type=str,help="Path to index to map against")
    parser.add_argument('-phred64',action="store_true",default=False,help = 'Use phred64 quality scores instead of phred33')
    parser.add_argument('-nomulti',action="store_false",default=True,help = 'Do not accept multimaps (bowtie option -a off)')
    parser.add_argument('-nofails',action="store_false",default=True,help = 'Do not log failed mappings (bowtie option --no-unal on)')
    parser.add_argument('-log',action="store_true",default=False,help = 'Create an explicit log of the mappings')
    parser.add_argument('-threads',type=int,default=4, help='[default = 4] Number of threads to use')
    parser.add_argument('-logname',type=str,default='batch_log.csv', help='[default = batch_log.csv] name of the log file')
    parser.add_argument('-suffix',type=str,default='mapped', help='[default = mapped] SAM file suffix')
    args = parser.parse_args()
    #
    if args.log == False:
        print ''
        print '\033[1;4;94mBatch Running Bowtie2 on current directory\033[0m'
        print ''
        start_time,count = time.asctime(),0
        for fyle in sorted(glob.glob('*.fastq')):
            print '\033[92mMapping: \033[0m{}'.format(fyle)
            bowtie_batcher(fyle, args.index, args.threads, args.suffix, args.phred64, args.nomulti,args.nofails,args.log)
            print ''
            count+=1
        print '\033[1;4;94mStart Time:\033[0m ',start_time
        print '\033[1;4;94mEnd Time:\033[0m ',time.asctime()
        print '\033[1;4;92mTotal <.fastq> mapped:\033[0m ',count
        print ''
    if args.log == True:
        bucket = []
        for fyle in sorted(glob.glob('*.fastq')):
           bucket.append(time.asctime())
           info,new_fyle,command = bowtie_batcher(fyle, args.index, args.threads, 
                                                  args.suffix, args.phred64, args.nomulti,args.nofails,args.log)
           
           bucket.extend([fyle,new_fyle,command,info,time.asctime()])
        write_log_file(bucket,args.logname)



if __name__ == '__main__': 
    main()

