#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import os
from optparse import OptionParser
import math
import glob
import threading
from itertools import islice


class subThread (threading.Thread):
    def __init__(self, threadID, function, *para):


        threading.Thread.__init__(self)
        self.para = para
        self.fun = function
        self.threadID = threadID
        self.name = "Thread-"+self.threadID
    def run(self):
        self.fun(*self.para)



#load reactivitiy file
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


#Read tab delimited file
def read_t_file(in_file):
    result = [];
    with open(in_file, 'r') as f:
        for aline in f.readlines():
            temp = [];
            tline = aline.strip();
            tl = tline.split('\t');
            for i in range(0, len(tl)):
                temp.append(tl[i].strip());
            result.append(temp);
    return result



#create temporary sequence file (fasta) for folding
def create_temp_seq(id_s, seqs, path):
    with open(os.path.join(path,id_s+"_temp.fa"), 'w') as fh:        
        fh.write('>'+id_s)
        fh.write('\n')
        if id_s in seqs:
            fh.write(seqs[id_s]+"\n")
            return True
        else:
            id_n.append(id_s)
            return False

#create temporary constraint file (fasta) for folding
def create_temp_react(id_s, react, path, shift, thres, mode):
    if id_s not in react:
        cons_n.append(id_s)
        os.system("rm "+os.path.join(path, id_s+"_temp.fa"))
        return False
    else:
        if not thres:
            with open(os.path.join(path, id_s+"_constraint.txt"), 'w') as fh:
                for j in range(0, len(react[id_s])-shift):
                    if react[id_s][j]!='NA':
                        fh.write(str(j+1))
                        fh.write('\t')
                        fh.write(str(react[id_s][j])+"\n")
        else:
            if mode == "1":
                with open(os.path.join(path, id_s+"_constraint.txt"), 'w') as fh:
                    ssn = []
                    for j in range(0, len(react[id_s])-shift):
                        if react[id_s][j]!='NA':
                            if float(react[id_s][j])>=thres:
                                ssn.append((j+1))
                    fh.write("DS:\n-1\nSS:\n")
                    for j in range(len(ssn)):
                        fh.write(str(ssn[j])+"\n")
                    fh.write("-1\nMod:\n-1\nPairs:-1 -1\nFMN:\n-1\nForbids:\n-1 -1\n")
            else:
                with open(os.path.join(path, id_s+"_temp.fa"), 'a') as fh:
                    for j in range(0, len(react[id_s])-shift):
                        if react[id_s][j]!='NA':
                            if float(react[id_s][j])>=thres:
                                fh.write('x')
                            else:
                                fh.write('.')
                        else:
                            fh.write('.')
                    for i in range(shift):
                        fh.write('.')
                                
                    
                    
        return True

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#Fold RNAs to predict structures
def fold_RNA(mode, p_type, id_s, temperature, slope, intercept, path, thres):
    t = id_s
    ct_path = os.path.join(path, "CT")
    ps_path = os.path.join(path, "PS")
    dbn_path = os.path.join(path, "DBN")
    bp_path = os.path.join(path, "Base-pairing_probability")

    
    if mode == '1':
        if not p_type:
            if not par:
                os.system(os.path.join(rp, "Fold")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" -md "+str(max_d)+" "+mf+os.path.join(ct_path, t+"_"+temperature+"_silico.ct")+" > /dev/null")
                os.system(os.path.join(rp, "draw")+" "+os.path.join(ct_path, t+"_"+temperature+"_silico.ct")+" "+os.path.join(ps_path, t+"_"+temperature+"_silico.ps")+" > /dev/null")
            else:
                os.system(os.path.join(rp, "partition")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" "+os.path.join(bp_path, t+"_"+temperature+"_silico.pfs")+" > /dev/null")
                os.system(os.path.join(rp, "ProbabilityPlot")+" "+os.path.join(bp_path, t+"_"+temperature+"_silico.pfs")+" -t "+os.path.join(bp_path, t+"_"+temperature+"_silico.bp")+" > /dev/null")
                os.system(os.path.join(rp, "ProbabilityPlot")+" "+os.path.join(bp_path, t+"_"+temperature+"_silico.pfs")+" "+os.path.join(ps_path, t+"_"+temperature+"_silico.ps")+" > /dev/null")
                ld = file_len(os.path.join(bp_path, t+"_"+temperature+"_silico.bp"))
                if ld > 2:
                    os.system(os.path.join(rp, "MaxExpect")+" "+os.path.join(bp_path, t+"_"+temperature+"_silico.pfs")+" "+os.path.join(ct_path, t+"_"+temperature+"_silico.ct")+" > /dev/null")

                os.system("rm "+os.path.join(bp_path, t+"_"+temperature+"_silico.pfs"))
        else:
            if not par:
                if not thres:
                    os.system(os.path.join(rp, "Fold")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" -md "+str(max_d)+" "+mf+"-sh "+os.path.join(path,id_s+"_constraint.txt")+" -si "+intercept+" -sm "+slope+" "+os.path.join(ct_path, t+"_"+temperature+"_restraint.ct")+" > /dev/null")
                    
                else:
                    os.system(os.path.join(rp, "Fold")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" -md "+str(max_d)+" "+mf+"-c "+os.path.join(path,id_s+"_constraint.txt")+" "+os.path.join(ct_path, t+"_"+temperature+"_restraint.ct")+" > /dev/null")
                os.system(os.path.join(rp, "draw")+" "+os.path.join(ct_path, t+"_"+temperature+"_restraint.ct")+" "+os.path.join(ps_path, t+"_"+temperature+"_restraint.ps")+" > /dev/null")
                    
        
            else:
                if not threshold:
                    os.system(os.path.join(rp, "partition")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" "+"-sh "+os.path.join(path,id_s+"_constraint.txt")+" -si "+intercept+" -sm "+slope+" "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs")+" > /dev/null")

                else:
                    os.system(os.path.join(rp, "partition")+" "+os.path.join(path,id_s+"_temp.fa")+" -T "+temperature+" "+"-c "+os.path.join(path,id_s+"_constraint.txt")+" "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs")+" > /dev/null")
                os.system(os.path.join(rp, "ProbabilityPlot")+" "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs")+" -t "+os.path.join(bp_path, t+"_"+temperature+"_restraint.bp")+" > /dev/null")
                os.system(os.path.join(rp, "ProbabilityPlot")+" "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs")+" "+os.path.join(ps_path, t+"_"+temperature+"_restraint.ps")+" > /dev/null")
                    
                ld = file_len(os.path.join(bp_path, t+"_"+temperature+"_restraint.bp"))
                if ld > 2:
                    os.system(os.path.join(rp, "MaxExpect")+" "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs")+" "+os.path.join(ct_path, t+"_"+temperature+"_restraint.ct")+" > /dev/null")
                    
                os.system("rm "+os.path.join(bp_path, t+"_"+temperature+"_restraint.pfs"))
                          
                          
    else:
        if not p_type:
            suffix = "silico"
        else:
            suffix = "restraint"
        if par:
            pr = " -p"
        else:
            pr = ""
        if not p_type:
            os.system(os.path.join(vp, "RNAfold")+pr+" < "+os.path.join(path,id_s+"_temp.fa")+" -T "+str(float(temperature)-273.15)+" > "+os.path.join(dbn_path, t+"_"+temperature+"_"+suffix+".dbn"))
        else:
            if thres:
                os.system(os.path.join(vp, "RNAfold")+pr+" < "+os.path.join(path,id_s+"_temp.fa")+" -T "+str(float(temperature)-273.15)+" > "+os.path.join(dbn_path, t+"_"+temperature+"_"+suffix+".dbn"))
            else:
                os.system(os.path.join(vp, "RNAfold")+pr+" < "+os.path.join(path,id_s+"_temp.fa")+" -T "+str(float(temperature)-273.15)+" --shape="+os.path.join(path,id_s+"_constraint.txt")+" > "+os.path.join(dbn_path, t+"_"+temperature+"_"+suffix+".dbn"))
            
                             
        if not par:
            os.system("mv "+id_s+"_ss.ps "+os.path.join(ps_path, t+"_"+temperature+"_"+suffix+".ps"))
        else:
            os.system("mv "+id_s+"_dp.ps "+os.path.join(ps_path, t+"_"+temperature+"_"+suffix+".ps"))
            os.system("rm "+id_s+"_ss.ps")
        
#Remove temporary files        
def remove_files(path, mode, id_s, p_type, thres):
    os.system("rm -f "+os.path.join(path,id_s+"_temp.fa"))
    if p_type:
        if mode == '1':
            os.system("rm -f "+os.path.join(path,id_s+"_constraint.txt"))
        else:
            if not thres:
                os.system("rm -f "+os.path.join(path,id_s+"_constraint.txt"))


def batch_predict(mode, predict_type, seqs, idss, temperature, slope, intercept, shift, threshold, react, output_directory, ml, ma):

    for i in range(len(idss)):
        idt = idss[i]
        flag = 0
        if len(seqs[idt])<ml:
            rts.append(idt)
        else:
            if len(seqs[idt])>ma:
                rtl.append(idt)
                
            else:
                if create_temp_seq(idt, seqs, output_directory):
                    flag = 1
                if predict_type and flag == 1:
                    if not create_temp_react(idt, react, output_directory, shift, threshold, mode):
                        flag = 0
                if flag == 1:       
                    fold_RNA(mode, predict_type, idt, temperature, str(slope), str(intercept), output_directory, threshold)
                    remove_files(output_directory, mode, idt, predict_type, threshold)
            
#separate the ID list into groups with equal distrubtion of length
def group_separate(id_sort, pro):
    r = []
    for i in range(pro):
        r.append([])
    for i in range(len(id_sort)):
        t = id_sort[i][0]
        g = i%pro
        r[g].append(t)
    return r
        

#progress bar
def update_meter_S(linecount,threshold):
    #print("dd")
    if linecount % threshold == 0:
        #sys.stdout.write('\033[1;94m'+u'\u25A3'+'\033[0m')
        sys.stdout.write('=')
        sys.stdout.flush()

def progress_bar(output_directory, thres):
    sys.stdout.write('[')
    sys.stdout.flush()
    ld = 0
    while(True):
        d = glob.glob(os.path.join(output_directory, 'PS', '*.ps'))
        lnew = len(d)+len(rts)+len(rtl)+len(id_n)+len(cons_n)
        if lnew > ld:
            update_meter_S(lnew, thres)
            ld = lnew
        if stop:
            sys.stdout.write(']')
            print("\n")
            break


def main():
    parser = argparse.ArgumentParser(description='Predict RNA structure from sequence (with restraints)\nAll output files are in the output_files_[temperature] folder', formatter_class=RawTextHelpFormatter)
    parser.add_argument('id_f', metavar="<RNA_IDs>",help='List of RNAs(IDs) to predict RNA structure')
    parser.add_argument('seq', metavar="<RNA_sequence>",help='File contains all the RNA sequences of the RNA IDs to be predicted (fasta format), e.g. Reference Transcriptome')
    parser.add_argument('-T', '--Temperature', dest="temperature", default = 310.15, help='The temperature under which the RNA structures are predicted [Default 310.15 K]', type = float)
    parser.add_argument('program', metavar="<Program_mode>",help='The program used for RNA structure prediction: 1 for RNAstructure; 2 for Vienna package')
    parser.add_argument("-r", "--reactivity", dest="Constraint_file", default = None, help='Reactivity file (.react file) used as constrants from prediction')
    parser.add_argument("-th", "--threshold", dest="thres", default = None, type = float, help='Threshold for reactivity. Any nucleotide with reactivity over threshold will be set as single stranded. Once threshold is set, reactivities will be converted into hard contraints')
    parser.add_argument("-sm", "--slope", dest="slope", default = 1.8, type = float, help='Slope for RNA structure prediction with restrains [Default 1.8] (Parameter for RNAstructure, only for prediction using RNAstructure)')
    parser.add_argument("-si", "--intercept", dest="intercept", default = -0.6, type = float, help='Slope for RNA structure prediction with restrains [Default -0.6](Parameter for RNAstructure, only for prediction using RNAstructure)')
    parser.add_argument("-sht", "--shift", dest="N", default = 0, type = int, help="Ignore the reactivities on the last N nucleotide of each RNA to be predicted [Default 0]")
    parser.add_argument("-minl", "--minumum_length", dest="MINL", default = 10, type = int, help="The minumum length of the RNA required for folding [Default 10]")
    parser.add_argument('result_s', nargs='?', metavar="<short_RNA_file>",help='Output IDs with corresponds to RNAs that are too short for folding [Default: Reads_too_short.txt]')
    
    parser.add_argument("-maxl", "--maximum_length", dest="MAXL", default = 5000, type = int, help="The maximum length of the RNA required for folding [Default 5000]")
    parser.add_argument('result_l', nargs='?', metavar="<long_RNA_file>",help='Output IDs with corresponds to RNAs that are too long for folding [Default: Reads_too_long.txt]')
    parser.add_argument("-par", "--PAR", action="store_true", help="Calculate partition function and output base-pairing probabilities instead of RNA structures")
    parser.add_argument("-mul", "--multiple", action="store_true", help="Output multiple predicted RNA structures instead of just outputing the MFE structure (only for RNAstructure prediction)")
    
    parser.add_argument("-p", "--process", dest="P", default = 1, type = int, help="Number of threads for parallel computing [Default 1]")
    parser.add_argument("-md", "--maxdistance", dest="MD", default = 99999, type = int, help="Specify a maximum pairing distance between nucleotides [Default: no restraint]")
    
    args = parser.parse_args()

    id_file = args.id_f
    seq_file = args.seq
    predict_type = args.Constraint_file
    temperature = args.temperature
    mode = args.program
    threshold = args.thres
    slope = args.slope
    intercept = args.intercept
    shift = args.N
    mul = args.multiple
    ml = args.MINL
    ma = args.MAXL
    mdd = args.MD
    
    pro = args.P

    global par
    par = args.PAR


    if args.result_s:
        result_file = args.result_s
    else:
        result_file = "RNAs_too_short.txt"

    if args.result_l:
        result_file_max = args.result_l
    else:
        result_file_max = "RNAs_too_long.txt"


    result_file_id_n = "ID_not_in_reference.txt"
    result_file_cons = "ID_not_in_reactivity.txt"

    global mf 
    if mul:
        mf = ""
    else:
        mf = "-mfe "

    global max_d
    max_d = mdd
    


    #shift = max(int(shift)-1, 0)
    program = {'1':"RNAstructure", '2':'Vienna_package'}

    if mode!= '1' and mode!='2':
        sys.stderr.write("Error: Mode selection error!\n")
        sys.exit()

    react = {}
    if predict_type:
        react = read_derived_reactivities(predict_type)

    if par:
        pn = "_partition_function"
    else:
        pn = ""
    
    syspath = os.getcwd()
    temperature = str(temperature)

    if mode == '2':
        mfe = ""
    
    if not predict_type:
        output_name = "silico_"+os.path.basename(id_file)+"_"+temperature+"_"+os.path.basename(seq_file)+"_"+program[mode]+mf.strip()+"_md_"+str(max_d)
    else:
        if threshold:
            output_name = os.path.basename(predict_type)+"_"+os.path.basename(id_file)+"_"+temperature+"_"+os.path.basename(seq_file)+"_"+program[mode]+mf.strip()+"_sht_"+str(shift)+pn+"_md_"+str(max_d)+"_thres_"+str(threshold)
        else:
            output_name = os.path.basename(predict_type)+"_"+os.path.basename(id_file)+"_"+temperature+"_"+os.path.basename(seq_file)+"_"+program[mode]+mf.strip()+"_sht_"+str(shift)+pn+"_md_"+str(max_d)
    output_directory = os.path.join(syspath, "output_files_"+output_name)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        sys.stderr.write("Error: Directory exists!\n")
        sys.exit()


    global rp
    global vp
    

    try:
        with open("set_path.txt", 'r') as f:
            for aline in f.readlines():
                line = aline.strip()
                p = line.split("=")[0].strip()
                if len(line.split("=")) > 1:
                    loc = line.split("=")[1].strip()
                else:
                    loc = ""
                if p == "RNAstructure":
                    rp = loc
                if p == "Vienna package":
                    vp = loc
    except:
        rp = ""
        vp = ""

    idss = read_t_file(id_file)
    idsu = set()
    for i in range(len(idss)):
        idsu.add(idss[i][0])
    sequences = SeqIO.parse(seq_file, 'fasta')
    seqs = {}
    leng = {}
    for seq in sequences:
        seqs[seq.id] = str(seq.seq)
        leng[seq.id] = len(str(seq.seq))

    ids_tup = []
    for tt in idsu:
        ids_tup.append((tt, leng[tt]))

    ids_sorted = sorted(ids_tup, key=lambda tup: tup[1])
    groups = group_separate(ids_sorted, pro)

    #print(ids_sorted)
    #print(groups)

    global rts
    rts = []

    global rtl
    rtl = []

    global id_n
    id_n = []

    global cons_n
    cons_n = []


    ps_path = os.path.join(output_directory, "PS")


    if not os.path.exists(ps_path):
        os.makedirs(ps_path)

    if par and mode == '1':
        bp_path = os.path.join(output_directory, "Base-pairing_probability")
        if not os.path.exists(bp_path):
            os.makedirs(bp_path)
                          
    if mode == '1':
        ct_path = os.path.join(output_directory, "CT")
        if not os.path.exists(ct_path):
            os.makedirs(ct_path)
    else:
        dbn_path = os.path.join(output_directory, "DBN")
        if not os.path.exists(dbn_path):
            os.makedirs(dbn_path)
        

    threads = []


    for i in range(len(groups)):
        idts = groups[i]
        #print(len(idts))
        thread = subThread(str(i), batch_predict, mode, predict_type, seqs, idts, temperature, slope, intercept, shift, threshold, react, output_directory, ml, ma)
        thread.start()
        threads.append(thread)

    i = i+1
    global stop
    stop = False
    thread = subThread(str(i), progress_bar, output_directory, max(1, len(list(idsu))/100))
    thread.start()

    


        
            
   

    for i in range(len(threads)):
        threads[i].join()

    stop = True

    with open(os.path.join(output_directory, result_file), 'w') as h:
        for i in range(len(rts)):
            h.write(rts[i]+"\n")


    with open(os.path.join(output_directory, result_file_max), 'w') as h:
        for i in range(len(rtl)):
            h.write(rtl[i]+"\n")

    if len(id_n)>0:
        with open(os.path.join(output_directory, "ID_not_in_reference.txt"), 'w') as h:
            for i in range(len(id_n)):
                h.write(id_n[i]+"\n")

    if len(cons_n)>0:
        with open(os.path.join(output_directory, "ID_not_in_reactivity.txt"), 'w') as h:
            for i in range(len(cons_n)):
                h.write(cons_n[i]+"\n")

        

    

if __name__ == "__main__" :
    main()
 
