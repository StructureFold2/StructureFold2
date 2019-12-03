# StructureFold2 <img src='assets/sf2_logo.png' align='right' width='400px' />
<br><br><br><br>

Please review the included manual, SF2_Manual.pdf, for information on using StructureFold2.

**Dependencies**
+ [Python 2.7.X](https://www.python.org/)
+ [BioPython](https://biopython.org/)
+ [Numpy](https://numpy.org/)
+ [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
+ [SAMtools](http://samtools.sourceforge.net/)
+ [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

**Recommended**
+ [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
+ [RNAStructure](https://rna.urmc.rochester.edu/RNAstructure.html)
+ [R](https://www.r-project.org/)

## Updates and Errata

**Planned Updates**
* <.ct> files generated under MFE settings now containing a DeltaG value, thus...
    * make_standed_csv.py and make_MFE_csv.py will be merged into one module.
    * make_PPV_csv.py will also be rolled into this module; structure_statistics.py.
    * Should be able to report/consolidate any number of directories of <.ct> files ito a single file.
