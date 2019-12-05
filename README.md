# StructureFold2 <img src='assets/sf2_logo.png' align='right' width='400px' />

Please review the included manual, SF2_Manual.pdf, for information on using StructureFold2.

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
structure_statisics.py has replaced make_standed_csv.py and make_MFE_csv.py, with the added
functionality of being able to combine multiple directories of <.ct> (connectivity tables) files into a single report. 
Older versions of RNAStructure did not include a DeltaG value for folds performed under MFE settings, but the more consistent
reporting of <.ct> files allows both strandedness and DeltaG to be extracted from such MFE folds, thus makine one module sufficient
for both tasks. 

#### Usage
```
Summarizes or compares MFE <.ct> files

optional arguments:
  -h, --help      show this help message and exit
  -d D [D ...]    CT directory/directories
  -mode {R,F}     Raw/Fused statistics
  -name NAME      Output file name
  -offset OFFSET  Number of Underscores in Transcript Names
```

## Planned Updates

* make_PPV_csv.py will be merged into structure_statisics.py
* coverage_overlap.py will be merged into rtsc_coverage.py
* Features that are exclusive to <.react> files have been requested to be available for <.rtsc>,
or vice-versa. I.E, rather than making a separate module for react_correlation to complement rtsc_correlation,
both features will be included in a single module; modules that work on both will be prefixed with 'rx'. So react_correlation and rtsc_correlation
would become rx_correlation, thus deprecating the old module(s) entirely.

