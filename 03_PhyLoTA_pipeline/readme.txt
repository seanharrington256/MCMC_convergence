Properties of Markov Chain Monte Carlo Performance Across Many Empirical Alignments

Acquisition and processing of phylogenetic datasets

This document describes the methods and scripts we used to obtain datasets from the phylota database http://sirloinpope.com/ (Sanderson et al. 2008) and conduct phylogenetic analyses. Briefly, fasta files were cleaned and renamed, reads were MUSCLE aligned, Gblocks was performed, jmodeltest was run and MrBayes blocks were written. These dataset files were then uploaded to a supercomputer cluster for analysis. Once complete, output files were downloaded for diagnostics and evaluation.

Preparing files for the pipeline
These scripts are located in the PhyLoTA_dataset_download directory.

Datasets were downloaded using a wget scripts to pull all datasets matching predefined criteria for version 184 get184datasets. In this case, all datasets where n taxa between 25 and 250, as determined using mySQL. Once raw datasets were downloaded, html headers were removed and unreadable characters by downstream programs were replaced by '_' using the script 001_cleanfiles.py. Finally, files were renamed to a simple common format encompassing the pertinent details (version number, root taxon ID, and cluster ID; 002_renamefiles. 27,307 total datasets were downloaded. Of these, 35 were empty (Zero bytes), and are removed, leaving 27,272 datasets for analysis.

General Pipeline Procedures
The following scripts are located within the PhyLoTA_analysis_pipeline directory.

Initial step - Writing directories
The bash script 01_directoryMaker writes necessary directories for the pipeline. Once fasta input files are copied in the 01_FASTA directory, the entire General Pipeline Procedures can be run in sequence by evoking the 00_masterScript file in the 00_SCRIPTS directory.

Pipeline Step 1 - Removing duplicate taxa
Fasta files were placed in the 01_FASTA folder and duplicate taxa were removed using a python script 02_duplicateRemover which reads the original fasta file, and writes a new file using only the first unique taxon ID's encountered. Thus, the new file will contain only unique taxa. These single taxon fasta files carry the .output extension.

Intermediate step - Moving output files
The bash script 03_outputMover moves the original fasta files (with duplicate taxa) to the /ORIGINAL_FASTA_FILES directory. The .output files (single-taxon files) are moved to the ../02_1TAXON_FASTA directory. Prior to moving files to the 02_1TAXON_FASTA directory, files may be manually inspected to ensure that each of the files contain >1 taxon. Otherwise, Gblocks will crash.

Pipeline Step 2 - Performing MUSCLE alignments
The bash script 04_MUSCLE_aligner performs MUSCLE alignments (Edgar 2004) on all files in the 02_1TAXON_FASTA directory with a .output extension. MUSCLE alignments are performed using default settings.

Intermediate step - Moving alignment files, and shortening taxa names
The bash script 05_alignmentMover moves all files which carry a .fasta extension (alignment files) to ../03_../MUSCLE_ALIGNMENTS. Files carrying the .output file extension are moved to ../02_1TAXON_FASTA/SINGLE_SPECIES_FASTA_FILES/.

The perl script 05.1_fasta_shorten.pl (Estill and Bennetzen 2009) changes the headers in the fasta files to give shorter names. This is necessary because Gblocks will crash if fasta headers are too long. In this case, all names are shortened to 50 characters by running the following command:
./05.1_fasta_shorten.pl -l 50 -i ../03_MUSCLE_ALIGNMENTS/ -o ../03_MUSCLE_ALIGNMENTS/SHORT_ID

Pipeline Step 3 - Gblocks
The program Gblocks (Castresana 2000) eliminates poorly aligned positions and divergent regions of nucleotide alignment sequence. These base positions may not be homologous or may be saturated by many substitutions so it is appropriate to them prior to phylogenetic analysis.
There are two steps in the Gblocks program: 1- writing a paths file, and 2- running the Gblocks program. The Gblocks program must read a paths file (directory path to files to be analyzed) to performs the analysis.  The bash script 06_makePathsFile writes file paths using a coreutils package: greadlink, which returns canonicalized name of a set of files—i.e., their absolute pathnames. This can be installed by running:brew install coreutils. Then, the Gblocks program can be called using the 07_Gblocks script, which pulls the paths from the file paths in the 00_SCRIPTS folder.

Pipeline Step 4 - Nexus file conversion
Since jmodeltest requires a nexus input file, it is necessary to convert the fasta files to nexus file format. This is accomplished by running the script 08_GblocksNexusConveter. This script runs Seqmagick, a utility that converts between several phylogenetic file types. In this case, we are using it convert from fasta to nexus.
Seqmagick can be installed by running pip install seqmagick.
Any fasta files in which Gblocks could not identify any blocks (i.e., empty sequences) are filtered out at this point, and not converted to nexus format.

Intermediate step - Moving nexus files
The bash script 09_nexusMover moves all nexus output files (.nex) to ../04_NEXUS_CONVERSION. It then moves all .fasta files (the alignment files) to ../03_MUSCLE_ALIGNMENTS/ALIGNMENTS. Finally, all Gblocks output files (.fas, .fas-e, and .htm) are moved to GBLOCKS_OUTPUT and the fas-e files are removed.

Pipeline Step 5 - Model testing
jModelTest (Guindon and Gascuel, 2003) is a tool that is used to carry out statistical selection of best-fit models of nucleotide substitution. The python script 10_jmodeltest.py implements this. Specifically, models are tested and evaluated using the following command line arguments:
-g 4 [include models with rate variation among sites, and set the number of categories to 4]
-i [include models with proportion invariable sites]
-f [include models with unequal base frequencies]
-t BIONJ [sets the base tree for likelihood calculations to Neighbor-Joining (NJ) topology for each model] (while NJ may be less accurate than a Maximum Likelihood [ML] tree topology search for each model, ML is computationally prohibitive for the amount of data considered here).
-s 3 [sets the number of substitution schemes to 3: JC/F81, K80/HKY, SYM/GTR, i.e., those implemented in MrBayes]
-AICc [select best-fit models using the corrected Akaike Information Criterion]

Pipeline Step 6 - Write Bayes block
The python script 11_mbb.py appends a Bayes block to the end of the nexus file. This contains important information for running phylogenetic analysis in MrBayes (Huelsenbeck and Ronquist, 2001), such as: specifying the substitution model to be used; priors; number of runs and chains; generations to run; diagnostic, print, and sample frequencies; as well as number of runs to be discarded as burn-in. MrBayes uses Markov Chain Monte Carlo to sample from the posterior probability distribution, and metropolis-coupling to speed up convergence where three "heated" chains are in parallel to the "cold" chain. The "heated" chains are simply sampling the same probability distribution, but one that has been "heated" to "melt" or "flatten" the landscape of peaks (optima) defined by the (posterior) probability distribution. Sampling switches between the cold and heated chains so the cold chain can escape local optima when sampling along the flattened peaks of probability space and more reliably identify the globally optimum parameters.

Intermediate steps - Moving analysis files, writing SLURM files, writing directories
The bash script 12_analysisFilesMover moves jModelTest output files to ../04_NEXUS_CONVERSION/jmodeltestOutputFiles and copies nexus (*.nex) output files to ../05_ANALYSIS_FILES/NEXUS_FILES.

References
Castresana, J. 2000. Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Mol Biol and Evol 17:540-5

Edgar, R. C. 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Res 32:1792-97.

Estill, J. C. and J. L. Bennetzen 2009. The DAWGPAWS pipeline for the annotation of genes and transposable elements in plant genomes. Plant Methods. 5:8

Guindon S. and O. Gascuel. 2003. A simple, fast and accurate method to estimate large phylogenies by maximum-likelihood. Systematic Biology 52:696-704.

Huelsenbeck, J. P. and F. Ronquist. 2001. MRBAYES: Bayesian inference of phylogeny. Bioinformatics 17:754-755.

Sanderson, M. J., D. Boss, D. Chen, K. A. Cranston, and A. Wehe. 2008. The PhyLoTA Browser: processing GenBank for molecular phylogenetics research. Syst. Biol. 57:335-346.
