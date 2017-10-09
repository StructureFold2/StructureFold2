#!/usr/bin/env python

'''
StructureFold2 batch running script for post-processing <.sam> alignment files before generating <.rtsc> files
This is considered an important part of the StructureFold2 pipleline and should be used even if you chose not to use bowtie2 as your
aligner or used the other batch scripts in previous steps
This script serves two primary purposes
1)Filter out reads in <.sam> files which are either unmapped or improperly mapped
2)Filter out reads in <.sam> files with more than x mismatches, or a mismatch on base position 1 (improper RT stop)
'''

#Imports
import glob
import subprocess
import argparse

#Functions
def filter_sam_file(infile,out_suffix='filtered',max_mismatch=3,allow_bp_first=False):
    '''docstring'''
    exceedmismatch,firstmismatch,exceedandfirst,header=0,0,0,0
    #outfile = infile.strip('samtools.sam')+out_suffix+'.sam'
    outfile = infile.replace('samtools.sam',out_suffix+'.sam')
    with open(infile) as f,open(outfile,'w') as g:
        for line in f:
            if line.startswith('@'):
                header+=1
            
            flag1,flag2,flag_nmmd,tl=0,0,0,line.strip().split('\t')
            
            for i in range(len(tl)):
                #Determine whether number of mismatches is over limitation
                if tl[i].strip().find("NM:i:")!=-1:
                    n = tl[i].strip().split(':')[-1].strip()
                    if int(n) > max_mismatch:
                        flag1 = 1
                    flag_nmmd = 1
                
                #Determine whether the read has first nt mismatch
                if tl[i].strip().find("MD:Z:")!=-1:
                    md = tl[i].strip().split(':')[-1].strip()
                    if (int(md[0])) == 0:
                        flag2 = 1
                    flag_nmmd = 1
            
            if flag_nmmd == 1:
                # If don't allow 1st nt mismatch
                if not allow_bp_first:
                    if flag1 == 1 and flag2 == 1:
                        exceedandfirst+=1
                    elif flag1 == 1 and flag2 == 0:
                        exceedmismatch+=1 #number of reads with mismatches is over limitation excluding those with 1st mismatch
                    elif flag1 == 0 and flag2 == 1:
                        firstmismatch+=1  ##number of reads with 1st nt mismatch excluding those with mismatches over limitation
                else:
                    #If allows 1st nt mismatch
                    if flag1 == 1:
                        exceedmismatch+=1
                
                if flag1 == 0:
                    if allow_bp_first:
                        g.write(line)
                    elif flag2 == 0:
                        g.write(line)

    return outfile,exceedandfirst,exceedmismatch,firstmismatch,header

def filter_sam_file_turbo(infile,out_suffix='filtered'):
    '''docstring'''
    #outfile = infile.strip('samtools.sam')+out_suffix+'.sam'
    outfile = infile.replace('samtools.sam',out_suffix+'.sam')
    with open(infile) as f,open(outfile,'w') as g:
        for line in f:
            flag,flag2,tl=0,0,line.strip().split('\t')
            for i in range(len(tl)):
                
                if tl[i].strip().find("NM:i:")!=-1:
                    n,flag2 = tl[i].strip().split(':')[-1].strip(),1
                    if int(n) > 3:
                        flag = 1
                
                if tl[i].strip().find("MD:Z:")!=-1:
                    md,flag2 = tl[i].strip().split(':')[-1].strip(),1
                    if (int(md[0])) == 0:
                        flag = 1
            if flag2 == 0:
                continue
            if flag == 0:
                g.write(line)
    return outfile

def file_line_count(afyle):
    '''Gets how many lines are in a file via subprocess'''
    return int(subprocess.Popen(['wc', '-l',afyle], stdout=subprocess.PIPE).communicate()[0].strip().split()[0])

def run_samtools_turbo(samfyle):
    '''Runs samtools with default settings quickly'''
    #outfile= samfyle.strip('.sam')+'_samtools.sam'
    outfile = samfyle.replace('.sam','_samtools.sam')
    subprocess.call(' '.join(['samtools','view','-h','-F','20','-S',samfyle,'>',outfile]),shell=True)
    return outfile

def run_samtools(samfyle,bitflag):
    '''Runs samtools with the specified combination of bitflags, returns a tuple with the file lengths before and after processing'''
    #outfile= samfyle.strip('.sam')+'_samtools.sam'
    outfile = samfyle.replace('.sam','_samtools.sam')
    before_count = file_line_count(samfyle)
    command = ' '.join(['samtools','view','-h','-F',str(bitflag),'-S',samfyle,'>',outfile])
    procced_info = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, 
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    empty_data = procced_info.stdout.read()
    after_count = file_line_count(outfile)
    return outfile,before_count,after_count

def write_log_file(info_bucket,outfile):
    '''Takes the info collected from stout and makes a nice consistent <.csv>'''
    sub_buckets = [info_bucket[x:x+12] for x in range(0,len(info_bucket),12)]
    x_header = ','.join(['in_sam','sam_lines','filtered_lines','sam_filter_flag','flag_options',
                         'mismatches_and_first_mismatch','mismatches','first_mismatch',
                         'header_lines','max_mismatch','out_sam','out_sam_lines'])
    with open(outfile,'w') as g:
        g.write(x_header+'\n')
        for record in sub_buckets:
            out_fields = ','.join(str(x) for x in record)
            g.write(out_fields+'\n')

#This may come in handy if we need to make this script aware of the version of SAMtools
def get_samtools_version():
    '''Gets the current version of samtools in the path'''
    procced_info = subprocess.Popen('samtools', shell=True, stdin=subprocess.PIPE, 
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True).stdout.read()
    info_lines = procced_info.split('\n')
    target_line = info_lines[2]
    return float('.'.join(target_line.split()[1].split('.')[0:2]))

def main():
    parser = argparse.ArgumentParser(description='\033[1;4;94mBatch filter all <.sam> files within a given directory. This will make them ready to be converted into <.rtsc>\033[0m')
    parser.add_argument('-turbo',action='store_true',default=False,help='All filter options ignored, default settings, no log')
    parser.add_argument('-sam',default=None, nargs='+', help='Specific files to operate on')
    parser.add_argument('-keep_all',action='store_true',default=False,help='Keep all intermediate files')
    parser.add_argument('-keep_reverse',action='store_true',default=False,help='Keep mappings from the reverse strand')
    parser.add_argument('-remove_secondary',action='store_true',default=False,help='Remove secondary alignments')
    parser.add_argument('-allow_bp1_mismatch',action='store_true',default=False,help='Accept mappings with first base mismatches',dest='firstmm')
    parser.add_argument('-logname',type=str,default='filter_log.csv', help='[default = filter_log.csv] Name of the log file')
    parser.add_argument('-suffix',type=str,default='filtered', help='[default = filtered] filtered <.sam> file suffix')
    parser.add_argument('-max_mismatch',type=int, default=3, help='[default = 3] Maximum allowed mismatches/indels',dest='max_mm')
    args = parser.parse_args()
    
    #Turbo Mode, no logs, no options.
    if args.turbo == True:
        
        #Make list of files to process
        fyle_lyst = sorted(glob.glob('*.sam')) if args.sam == None else sorted(args.sam)
        for fyle in fyle_lyst:
            out1 = run_samtools_turbo(fyle)
            out2 = filter_sam_file_turbo(out1,args.suffix)
            if args.keep_all == False:
                subprocess.call(' '.join(['rm',out1]),shell=True)

    #Detailed Mode, slower, explicit log.
    else:
        #Assemble a bitflag to use for all samtools processing.
        unmapped = 4 #This currently must be set to 4.
        reverse = 16 if args.keep_reverse == False else 0
        secondary = 256 if args.remove_secondary == True else 0
        descriptions = {4:'unmapped',16:'reverse_strand',256:'secondary_alignments'}
        bitflag = str(sum([unmapped,reverse,secondary]))
        options_column = ' '.join([descriptions[flag] for flag in [x for x in [unmapped,reverse,secondary] if x!=0]])
        
        #Make a list of files to process, and a list to fill with information
        bucket = []
        fyle_lyst = sorted(glob.glob('*.sam')) if args.sam == None else sorted(args.sam)
        #Repeat both functions on the files
        for fyle in fyle_lyst:
            
            fyle_1,count_1,count_2 = run_samtools(fyle,bitflag)
            bucket.extend([fyle,count_1,count_2,bitflag,options_column])
            
            fyle_2,firstrsnp,snp,first,headlines= filter_sam_file(fyle_1,args.suffix,args.max_mm,args.firstmm)
            bucket.extend([firstrsnp,snp,first,headlines,str(args.max_mm),fyle_2,str(file_line_count(fyle_2))])
            
            if args.keep_all == False:
                subprocess.call(' '.join(['rm',fyle_1]),shell=True)
        
        write_log_file(bucket,args.logname)

if __name__ == '__main__':
    main()
