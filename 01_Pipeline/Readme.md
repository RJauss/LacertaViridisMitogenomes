## Required programs to run all scripts
- Trimmomatic / version 0.39
- MITObim / version 1.9.1
- Bowtie2 / version 2.3.5.1
- Samtools / version 1.10
- MAFFT / version 7.428
- Ruby Interpreter / version 2.6.0
- Ruby gems:
  - bio (BioRuby) / version 2.0.0
  - descriptive_statistics (Descriptive Statistics) / version 2.5.1
- RAxML / version 8.2.12
- MrBayes / version 3.2.6
- R (especially RScript)
  - R library APE

Remarks:
- most of the required programs must within search path (to start it from anywhere)
- some programs require to be installed at /bio-bin/ root path (Trimmomatic etc.) but scripts can be adjusted as usual.
- Ruby scripts `*.rb` have been written to better handle text based files and automate the processing where shell is not enough

## Tools required / used for manual steps
- Tablet v1.10.05.28
- MITOS http://mitos.bioinf.uni-leipzig.de/index.py
- ProtTest3 v3.4.2
- MEGA version X (tree visualization, can be done by any other prefered)

## Zip Package Layout

### Folder: Daten
- contains reference sequences
- contains nextgen raw data files (not yet included, too big: 5,68 GB)

### Folder: MITOS__fasta_for_mitos
- contains all final contig fasta files send to MITOS to anotate all gens

### Folder: MITOS_2019_08_09
- contains all results obtained from MITOS for each contig sent

### Folder: MITOS_2019_08_09_genes
- contains extracted gens required for each contig
- contains unaligned and aligned (MAFFT) alignments per gen
- contains generated ramxl complete fasta alignment with raxml partition file
- contains generated mrbayer complete nexus file with mrbayes partition block file
- contains the reconstruction sub folder containing the results from raxml and mrbayes

### Folder: NEXTGEN
- contains the `pipeline.sh` shell script to produce finally the contig files from nexgen data
- will create sub folders

### Folder: NEXTGEN / trimmed
- contains results of Trimmomatic with several trim sizes

### Folder: NEXTGEN / mitobim
- contains results of MITObim

### Folder: NEXTGEN / bowtie
- contains reference index db taken from reference sequence AM176577_1
- contains results of Bowtie2 and Samtools

### Folder: NEXTGEN / alignment
- contains fasta alignments (unaligned / aligned with MAFFT) with sequences for each trim size

### Scripts in Main Folder
- `fasta2nex.r` will be used to convert fasta files in nexus format as command line tool
- `gen_split.rb` will be used to split up the required gens from MITOS results
- `execute.rb` will be run after `gen_split.rb` to produce:
  - RAxML input files (fasta and partition)
  - MrBayes input files (nexus and partition)
  - automatic run and monitor RAxML and MrBayes execution
  - convert the results from RAxML and MrBayes to newick format (*.nwk)
- `dloop_stats.rb` checks all sequences in terms of motifs in dloop and reports results.