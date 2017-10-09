#!/usr/bin/env python
#StructureFold2 batch running script for using cutadapt on <.fastq> before aligning to a reference
#This is an accessory script and is not a core part of StructureFold2, but you may find it useful.

#Imports
import glob
import subprocess
import time
import argparse
import sys

def trim_fastq(afile,fiveprime,threeprime,minlen,minqual,maxlen,log_flag,suffix,nextseqflag):
    '''Issues the shell commands for each step in sequence'''
    if log_flag == False:
        out1,out2,out3 = [afile.split('.')[0]+'_'+str(x)+'.fastq' for x in ['1','2',suffix]]
        subprocess.call(' '.join(['cutadapt','-n','2','-g',fiveprime,'-o',out1,afile]),shell=True)
        subprocess.call(' '.join(['cutadapt','-a',threeprime,'-o',out2,out1]),shell=True)
        Q = ['-q',minqual] if nextseqflag == False else ['--nextseq-trim='+str(minqual)]
        M = ['-M',maxlen] if maxlen else []
        command = ['cutadapt','-m',minlen]+Q+M+['-o',out3,out2]
        #subprocess.call(' '.join(['cutadapt','-m',minlen,'-q',minqual,'-M',maxlen,'-o',out3,out2]),shell=True)
        subprocess.call(' '.join(command),shell=True)
        subprocess.call(' '.join(['rm',out1]),shell=True)
        subprocess.call(' '.join(['rm',out2]),shell=True)
    if log_flag == True:
        out1,out2,out3 = [afile.split('.')[0]+'_'+str(x)+'.fastq' for x in ['1','2',suffix]]
        FP_command = ' '.join(['cutadapt','-n','2','-g',fiveprime,'-o',out1,afile])
        TP_command = ' '.join(['cutadapt','-a',threeprime,'-o',out2,out1])
        Q = ['-q',minqual] if nextseqflag == False else ['--nextseq-trim='+str(minqual)]
        M = ['-M',maxlen] if maxlen else []
        TR_command = ' '.join(['cutadapt','-m',minlen]+Q+M+['-o',out3,out2])
        #TR_command = ' '.join(['cutadapt','-m',minlen,'-q',minqual,'-M',maxlen,'-o',out3,out2])
        FP_log = subprocess.Popen(FP_command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,close_fds=True).stdout.read()
        TP_log = subprocess.Popen(TP_command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,close_fds=True).stdout.read()
        TR_log = subprocess.Popen(TR_command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,close_fds=True).stdout.read()
        subprocess.call(' '.join(['rm',out1]),shell=True)
        subprocess.call(' '.join(['rm',out2]),shell=True)
        return FP_log+TP_log+TR_log

def dump_bucket(alist,outfyle):
    '''Dumps the logs of stdout to a file for concise record keeping.'''
    with open(outfyle,'w') as g:
        for line in alist:
            g.write(line)

def main():
    parser = argparse.ArgumentParser(description='Batch run cutadapt with given settings on all <.fastq> in a directory')
    parser.add_argument('-log',action="store_true",default=False,help = 'Create an explicit log of the trimming')
    parser.add_argument('-nextseq',action="store_true",default=False,help = 'Use NextSeq/NovaSeq quality scores')
    parser.add_argument('-fp',type=str, default='TGAACAGCGACTAGGCTCTTCA', help='[default = TGAACAGCGACTAGGCTCTTCA] 5\' adapter',dest='fp_adapt')
    parser.add_argument('-tp',type=str, default='GATCGGAAGAGCACACGTCTG', help='[default = GATCGGAAGAGCACACGTCTG] 3\' adapter',dest='tp_adapt')
    parser.add_argument('-minlen',type=str, default='20', help='[default = 20] minimum accepted sequence length', dest='min_len')
    parser.add_argument('-minqual',type=str, default='30', help='[default = 30] minumum accepted base quality',dest='min_qual')
    parser.add_argument('-maxlen',type=str, default=None, help='[default = None] maximum seq length',dest='max_len')#Make this optional
    parser.add_argument('-suffix',type=str,default='trimmed', help='[default = trimmed] trimmed <.fastq> file suffix',dest='suffix')
    parser.add_argument('-logname',type=str,default='trim_log', help='[default = trim_log] Name of the log file')
    args = parser.parse_args()
    #
    if args.log == False:
        print ''
        print '\033[1;4;94mStructure Fold2:\033[0;0;92m trim_fastqs.py\033[0m'
        print ''
        print '\033[94mMass trimming ALL <.fastq> files\033[0m'
        print ''
        start_time = time.asctime()
        fastqs,count = sorted(glob.glob('*.fastq')),0
        for fyle in fastqs:
            trim_fastq(fyle,args.fp_adapt,args.tp_adapt,args.min_len,args.min_qual,args.max_len,args.log,args.suffix,args.nextseq)
            count+=1
        print '\033[1;4;94mStart:\033[0m',start_time
        print '\033[1;4;94mFinished:\033[0m',time.asctime()
        print '\033[1;4;92mFASTQ Processed:\033[0m',count
        print ''
    else:
        bucket = []
        fastqs,count = sorted(glob.glob('*.fastq')),0
        for fyle in fastqs:
            start_time = time.asctime()+'\n'
            all_lines = trim_fastq(fyle,args.fp_adapt,args.tp_adapt,args.min_len,args.min_qual,args.max_len,args.log,args.suffix,args.nextseq)
            end_time = time.asctime()+'\n'
            bucket.append(start_time+all_lines+end_time)
        dump_bucket(bucket,args.logname+'.txt')

if __name__ == '__main__':
    main()

