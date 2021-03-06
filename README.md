# StructureFold2 <img src='assets/sf2_logo.png' align='right' width='400px' />

For information on using StructureFold2, Please consult the included manual, SF2_Manual.pdf.
If you use StructureFold2, please cite [_Tack et al._, 2018](https://www.sciencedirect.com/science/article/pii/S1046202317303535).


**Software Dependencies**
+ [Python 2.7.X](https://www.python.org/)
+ [BioPython](https://biopython.org/)
+ [Numpy](https://numpy.org/)
+ [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
+ [SAMtools](http://samtools.sourceforge.net/)
+ [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

**Recommended Software**
+ [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
+ [RNAStructure](https://rna.urmc.rochester.edu/RNAstructure.html)
+ [R](https://www.r-project.org/)

## Updates and Errata

### structure_statistics.py
Structure_statisics.py replaces make_standed_csv.py, make_MFE_csv.py, and make_PPV_csv.py with the added
functionality of being able to combine summaries or comparisons of multiple directories of <.ct> files 
(connectivity tables) generated by batch_fold_rna.py into a single report. 
Older versions of RNAStructure did not include a DeltaG value for folds performed under MFE settings, 
but the new more consistent reporting of <.ct> files allows both strandedness and DeltaG to be extracted from such 
MFE folds, thus making one module sufficient for all tasks involving <.ct> files. 

Fused mode (-mode F) summarizes one or multiple directories, generating a single <.csv> organized by 
transcript detailing transcript strandedness and DeltaG. When using multiple directories, the names of the directories 
must be different. If your transcript names contain underscores, input the number they contain (-offset). 
Raw mode (-mode R) instead operates on a single directory and generates a similar single <.csv>, instead organized by 
the files contained in the directory, and may be of use when it is hard to automatically parse transcript names. 
By default, 'NA' will be logged as the DetlaG for <.ct> files that do not contain a DeltaG value; this value is 
configurable (-na). PPV mode (-P) does a full directory-wise PPV calculation on every shared transcript between every 
two possible pairs of directories (-d) entered using [RNAStructure](https://rna.urmc.rochester.edu/RNAstructure.html)'s 
scorer program, which must be installed for this mode to work. All results will be summed into a single <.csv> report file.
In both Fused and PPV mode, the directories used do not have to contain perfectly symmetrical sets of <.ct> files.

**Usage**
```
Summarizes or compares MFE <.ct> files

optional arguments:
  -h, --help      show this help message and exit
  -d D [D ...]    CT directory/directories
  -mode {F,R,P}   Fused/Raw/PPV statistics
  -name NAME      Output file name
  -na NA          [default = NA] Null deltaG value
  -offset OFFSET  Number of Underscores in Transcript Names
```

### rtsc_coverage.py
The functionality of coverage_overlap.py has been combined into rtsc_coverage.py for convenience and
overall organization. Coverage overlap files are now generated concurently to the coverage calculation (-ol). 
Input file selection is now more precise for complex experimental designs (-f). In cases where
only a single condition is being probed, an overlap file is still used for restricting downstream
analyses to those transcripts meeting or exceeding the coverage threshold (-ot).

**Usage**
```
Creates a <.csv> of stop coverages from <.rtsc> files

positional arguments:
  index         <.fasta> file used to generate <.rtsc>

optional arguments:
  -h, --help    show this help message and exit
  -f F [F ...]  <.rtsc> files to operate on
  -bases BASES  [default = AC] Coverage Specificity
  -name NAME    Output file name
  -ol           Create an overlap file
  -ot OT        [default = 1.0] Overlap file threshold
  -on ON        Overlap file name
```

### rtsc_downscale.py
This module offers several avenues to downscale <.rtsc> files, creating
new <.rtsc> files with only a subset of all the reverse transcriptase stops. To
multiply the total reverse transcriptase counts of each base by a given ratio (-ratio), 
use FRACTIONAL mode; fractional stops do not interfere with the reactivity calculation.
Pseudo-random removal of stops can be performed using RANDOMREAD mode; each RT stop, regardless
of the base it occurs on has specificed chance (-ratio) to be retained. RANDOMPOSITION mode
continually selects a random position along each transcript until it can remove a stop, 
removing stops in this manner until a specified ratio (-ratio) of the total original 
number of stops remain on that transcript.

**Usage**
```
Downscales <.rtsc> files.

positional arguments:
  {FRACTIONAL,RANDOMREAD,RANDOMPOSITION}

optional arguments:
  -h, --help            show this help message and exit
  -f F [F ...]          Specific <.rtsc> to operate on
  -ratio RATIO          [default = 0.50] Fraction of RT stops to retain
  -restrict RESTRICT    Limit downscaling to these specific transcripts <.txt>
  -sort                 Sort output by transcript name
```

### react_motif.py
This module searches for instances of given motif(s) within two <.react> files,
thereby cataloging all of the reactivity change in and around those motifs. Both the 
control and experimental <.react> files must have been generated against the same <.fasta>, 
which also must be provided. A single motif can be entered on the command line (i.e. TATTA),
or a flat text file containing multiple motifs, one per line, can be entered. Bases directly
up or down stream of the motif can also be included in the output (-fp, -tp, default 5). 
In addition to the default <.csv> output, a new <.fasta> and corresponding <.react>s of the 
motif containing regions can be created (-fastaout, -reactout), where these files could be 
useful for guiding folding algorithms to regions containing motifs of interest. This module 
replaces react_static_motif.py entirely; all perumations of wildcard containing motifs
are now logged to a single file (<.csv>, <.fasta>,<.react>), and these files are written to 
a new directory, thereby enhancing organization.

**Usage**
```
Searches for reactivity differences around given motifs

positional arguments:
  control             control <.react> file
  experimental        experimental <.react> file
  fasta               <.fasta> to pull sequences from
  motif               Input file or motif

optional arguments:
  -h, --help          show this help message and exit
  -fp FP              [default = 5] Bases to include 5' of the motif
  -tp TP              [default = 5] Bases to include 3' of the motif
  -restrict RESTRICT  <.txt > Limit analysis to these specific transcripts
  -outdir OUTDIR      [default = motif_out] Out Directory
  -fastaout           Write windows to <.fasta> format as well
  -reactout           Write accompanying <.react> files as well
```

### react_heat_correct.py
Probing reagents can be more reactive at higher temperatures, thus one may wish to 
normalize this effect out from any two pairwise <.react>s derived from different temperatures. 
For any two such <.react> files generated against the same transcriptome, the sum of all base 
reactivties will be summed for both files, and used to scale the the lower/higher temperature file
up/down respectively, such that the total amount of reactivity on all bases in both files is the same. 

**Usage**
```
Corrects two <.react>s for differential temperature

positional arguments:
  lower           lower temp <.react> file
  higher          higher temp <.react> file

optional arguments:
  -h, --help      show this help message and exit
  -suffix SUFFIX  [default = corrected] Suffix for out files
```

### react_composition.py
Excessive amounts of wildcards used in react_motif.py result in a very computationally 
expensive and slow search, when the targeted motif may actually not be specific. If the 
user just wants any region of length n with a given base composition, for example all 
stretches of only G and C of over 10 bases in length, a different search method may be 
better suited for the task. react_composition.py addresses this demand, and functions in a similar fashion 
to react_motif.py. Instead of providing a motif, or a list of motifs with wildcards, windows
of a given length (-size) containing at least a specified percentage (-perc) of given bases 
(-bases) are found, and all reactivity changes in and around them are resolved and subsequently 
logged to a <.csv>. As with react_motif.py, these windows may be logged to fasta (-fastaout) format
and/or as react files.

**Usage**
```
Finds reactivity differences in windows of a given size and composition

positional arguments:
  control             control <.react> file
  experimental        experimental <.react> file
  fasta               <.fasta> to pull sequences from
  bases               Bases to query (i.e.'GC' or 'AT')

optional arguments:
  -h, --help          show this help message and exit
  -size SIZE          [default = 8] Size of window
  -perc PERC          [default = 1.0] Percent specified bases
  -unique             [default = True] Remove overlapping windows
  -fp FP              [default = 5] Bases to include 5' of the motif
  -tp TP              [default = 5] Bases to include 3' of the motif
  -restrict RESTRICT  <.txt > Limit analysis to these specific transcripts
  -fastaout           Write windows to <.fasta> format as well
  -reactout           Write accompanying <.react> files as well
```

### react_multi_motif.py
Instead of searching for areas of interest by composition (react_composition) or a single specified motif (react_motif), react_multi_motif
searches for areas featuring any of several (-mn, default = 3) smaller non-overlapping sub-motifs within a given window (-mw, default=50), 
otherwise sharing many of the settings and features similar modules. The -unique setting does not check if the 5'/3' regions 
flanking these multi-motifs overlaps other meta-motifs or their respective flanking regions, rather it only prevents 
the core multi-motifs themselves from overlapping.

**Usage**
```
Searches for reactivity differences around given multi-motifs

positional arguments:
  control             control <.react> file
  experimental        experimental <.react> file
  fasta               <.fasta> to pull sequences from
  motifs              multi-motif components

optional arguments:
  -h, --help          show this help message and exit
  -unique             [default = True] Remove overlapping multi-motifs
  -mn MN              [default = 3] Number of multi-motifs required
  -mw MW              [default = 50] Query window size for multi-motifs
  -FP FP              [default = 30] Bases 5' of multi-motifs
  -TP TP              [default = 30] Bases 3' of multi-motifs
  -restrict RESTRICT  <.txt > Limit analysis to these specific transcripts
  -fastaout           Write windows to <.fasta> format as well
  -reactout           Write accompanying <.react> files as well
```

## Planned Updates
* Features that are exclusive to <.react> files have been requested to be available for <.rtsc>,
or vice-versa. I.E, rather than making a separate module for react_correlation to complement rtsc_correlation,
both features will be included in a single module; modules that work on both will be prefixed with 'rx'. 
So react_correlation and rtsc_correlationwould become rx_correlation, thus 
deprecating the old module(s) entirely.<br><br>
* Add support for [STAR](https://github.com/alexdobin/STAR) into fastq_mapper.py. STAR is **much** faster than 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and will be integrated
into the SF pipeline in the future. STAR contains much better internal logging as well. 
Bowtie2 support will be maintianed for users that do not have access to a machine
with enough RAM run STAR.<br><br>
* Hardware guides (i.e. Linux workstation builds) for those labs
looking to get a machine to do bioinformatics may be created.<br><br>
* batch_fold_rna.py will be completely reworked to be more intuative and efficient. Alloting one thread per fold
works fine for a pool of smaller RNAs, but comitting only one thread on on larger RNAs really slows things down.<br><br>
* All motif modules may be merged into a single motif searching module with three modes (motif,composition,multi-motif)