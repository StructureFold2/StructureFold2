#Extra Tools for Reactivity Comparisons
import re
import itertools

class MotifEntry(object):
    '''Just holds all the components of an entry in a clean way'''
    def __init__(self,transcript,query,coords,*categories):
        self.transcript = transcript
        self.query = query
        self.pycoords = coords
        self.start = coords[0]+1
        self.end = coords[1]+1
        for sub_dict in categories:
            for k in sub_dict:
                setattr(self, k, sub_dict[k])

    def generate_key(self):
        return (self.transcript,self.query,self.start,self.end)

    def generate_seq(self):
        return self.fp_seq+self.seq+self.tp_seq

    def A_react(self):
        return self.fp_A+self.A+self.tp_A
 
    def B_react(self):
        return self.fp_B+self.B+self.tp_B

    def D_react(self):
        return self.fp_delta+self.delta+self.tp_delta

    def generate_name(self):
        return '_'.join([self.transcript,self.query,str(self.start),str(self.end)])

class MotifReport(object):
    '''Generates and Contains a full Motif Report'''
    def __init__(self,motif,fasta_dict,react_A,react_B,FP,TP,fill='-'):
        self.motif = motif
        self.queries = permute_motif(motif)
        self.FP = FP
        self.TP = TP
        self.size = FP+len(motif)+TP
        self.records = {}
        for query in self.queries:
            for transcript,sequence in fasta_dict.items():
                try:
                    A,B = react_A[transcript],react_B[transcript]
                    delta = subtract_reacts(B,A)
                    matches = generate_motif_coords_nested(sequence,query)
                    vectors = [A,B,delta,sequence]
                    for m in matches:
                        coreA,coreB,coreD,coreS = [core[m[0]:m[1]] for core in vectors]
                        fpA,fpB,fpD,fpS = [fp[max(m[0]-FP,0):m[0]] for fp in vectors]
                        tpA,tpB,tpD,tpS = [tp[m[1]:min(len(tp)-1,m[1]+TP)] for tp in vectors]
                        
                        if len(fpS) < FP:
                            fpfill = (FP-len(fpS))*fill
                            fpS = fpfill+fpS
                            fpA,fpB,fpD = [list(fpfill)+orange for orange in [fpA,fpB,fpD]]

                        if len(tpS) < TP:
                            tpfill = (TP-len(tpS))*fill
                            tpS = tpS+tpfill
                            tpA,tpB,tpD = [lemon + list(tpfill) for lemon in [tpA,tpB,tpD]]
                        
                        q_items = {'A':coreA,'B':coreB,'delta':coreD,'seq':coreS}
                        fp_items = {'fp_A':fpA,'fp_B':fpB,'fp_delta':fpD,'fp_seq':fpS}
                        tp_items = {'tp_A':tpA,'tp_B':tpB,'tp_delta':tpD,'tp_seq':tpS}
                        q_items.update(delta_metrics(coreD,'query'))
                        fp_items.update(delta_metrics(fpD,'fp'))
                        tp_items.update(delta_metrics(tpD,'tp'))
                        entry = MotifEntry(transcript,query,m,q_items,fp_items,tp_items)
                        self.records[entry.generate_key()] = entry
                except KeyError:
                    continue

def permute_motif(motif):
    '''Returns a list of all possible permutations of a DNA motif'''
    if all([base in dna for base in motif]):
            return [motif]
    else:
        motifs = []
        nightmare = [wildcards[base] for base in motif if base in wildcards]
        replacements = [iter(q) for q in itertools.product(*nightmare)]
        for replacement in replacements:
            permute = ''.join([base if base in dna else replacement.next() for base in motif])
            motifs.append(permute)
        return sorted(motifs)

def subtract_reacts(react_A,react_B):
    '''Subtracts B from A, i.e. A-B'''
    return [a-b if all([isinstance(a,float),isinstance(b,float)]) else 'NA' for a,b in zip(react_A,react_B)]

def add_reacts(react_A,react_B):
    '''Adds A and B, i.e. A+B'''
    return [a+b if all([isinstance(a,float),isinstance(b,float)]) else 'NA' for a,b in zip(react_A,react_B)]

def generate_motif_coords_basic(sequence,motif):
    '''Returns all non-overlapping coordinates of matching patterns'''
    return [match.span() for match in re.finditer(motif, sequence)]

def generate_motif_coords_nested(sequence,motif):
    '''Returns all coordinates of matching patterns'''
    return [match.span(1) for match in re.finditer('(?=('+motif+'))', sequence)]

def delta_metrics(delta_list,prefix):
    '''shiroyama'''
    i_numbs = filter(lambda i: isinstance(i,float),delta_list)
    p_numbs = filter(lambda p: p > 0,i_numbs)
    n_numbs = filter(lambda n: n < 0,i_numbs)
    a_numbs = [abs(a) for a in i_numbs]
    sums = [sum(l) for l in [i_numbs,p_numbs,n_numbs,a_numbs]]
    labels = ['change','increases','decreases','abs_change']
    prefixed_lables = ['_'.join([prefix,lab]) for lab in labels]
    return {numb:label for numb,label in zip(prefixed_lables,sums)}

#The standard DNA wildcard alphabet
wildcards = {'N':['A','C','G','T'],
    'R':['A','G'],'Y':['C','T'],
    'W':['A','T'],'S':['G','C'],
    'M':['A','C'],'K':['G','T'],
    'B':['G','C','T'],
    'H':['A','C','T'],
    'D':['A','G','T'],
    'V':['A','G','C']}

#Standard DNA and RNA bases
dna = ['A','C','G','T']
rna = ['A','C','G','U']
