# multiPhATE v.1.50
/MultiPhATE/ - multiPhATE2 (beta version) - This code is under development - Testing is appreciated!
For a production version of multiPhATE, use the code at https://github.com/carolzhou/multiPhATE.git.

This code was developed by Carol L. Ecale Zhou and Jeffrey Kimbrel at Lawrence Livermore National Laboratory.

THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD-3.pdf FOR DETAILS.

#### WHAT'S NEW?

1) multiPhATE inputs an optional user-provided (custom) gene-call set.
2) In comparing gene calls from multiple gene callers, multiPhATE now generates several new gene-call sets: superset, consensus, commoncore.
3) The new gene-call sets can be forwarded to the functional annotation processing.
4) The PhATE annotation pipeline now runs phmmer and hmmscan, in addition to jackhmmer.
5) PhATE now captures the raw hmm search output (ie, alignments).
6) multiPhATE now runs CompareGeneProfiles, a code that identifies gene similarities among genomes.
7) A genomics module computes gene and protein homology groups among all genomes input to the pipeline

#### ABOUT THE MULTI-PHATE PIPELINE DRIVER

MultiPhATE is a command-line program that runs the PhATE annotation code over user-specified phage genomes, then performs gene-by-gene comparisons among the genomes. The multiPhate.py code takes a single argument consisting of a configuration file (hereafter referred to as, "multiPhate.config"; use sample.multiPhate.config as starting point) and uses it to specify annotation parameters. Then, multiPhate.py invokes the PhATE pipeline for each genome. See below for the types of annotations that PhATE performs. If two or more genomes are specified by the user, then multiPhATE will run the CompareGeneProfiles code to identify corresponding genes among the genomes.

#### ABOUT THE PHATE PIPELINE

PhATE is a fully automated computational pipeline for identifying and annotating phage genes in genome sequence. PhATE is written in Python 3.7, and runs on Linux and Mac operating systems. Code execution is controled by a configuration file, which can be tailored to run specific gene finders and to blast sequences against specific phage- and virus-centric data sets, in addition to more generic (genome, protein) data sets. See below for the specific databases that are accommodated. PhATE runs at least one gene finding algorithm, then annotates the genome, gene, and protein sequences using nucleotide and protein blast flavors and a set of fasta sequence databases, and uses hmm searches (phmmer, jackhmmer) against these same fasta databases. It also run hmmscan against the pVOG profile database. If more than one gene finder is run, PhATE will provide a side-by-side comparison of the genes called by each gene caller. The user specifies the preferred gene caller, and the genes and proteins predicted by that caller are annotated using blast against the supporting databases (or, the user may specify one of the comparison gene sets: superset, consensus, or commoncore, for functional annotation). Classification of each protein sequence into a pVOG group is followed by generation of an alignment-ready fasta file. By convention, genome sequence files end with extension, ".fasta"; gene nucleotide fasta files end with, ".fnt", and cds amino-acid fasta files end with, ".faa".

#### ABOUT COMPARE-GENE-PROFILES and the GENOMICS MODULE

CompareGeneProfiles performs binary blast (NxN) of the genes from each genome against the genes from every other genome provided by the user. The code then identifies for each gene its mutual and non-mutual best hits against corresponding genes from each of the other genomes, and reports if no corresponding hit is found. For each binary genome-to-genome comparison, hits are ordered with respect to the query (reference, or first) genome. The Genomics module inputs the binary blast files from CompareGeneProfiles and computes genes and proteins that correspond across all the input genomes with respect to the reference genome. Ultimately, homology groups comprising each reference gene (or protein) and its corresponding genes, plus its homologs and their corresponding genes. Homology groups are output as fasta files and annotation files.

#### HOW TO SET UP MULTI-PHATE ON YOUR LOCAL MACHINE

First, create a working directory on your computer for running multiPhATE. Then, acquire the multiPhATE package from github. This can be done either by downloading a zip file directly from the multiPhATE repository, or by cloning the repository. The first method is recommended, but the second is certainly an option:

=) To download the zip file:  Use a browser and navigate to https://github.com/carolzhou/multiPhATE2. Press the green button "Clone or download", and download the zip file. Then, unzip the package in your working (main execution "multiPhate") directory.

```
$ cd myMultiphateDir

$ unzip multiPhate2-master.zip
```

=) To clone from github:  Acquire git from https://git-scm.com/downloads. Naviate to your working (main execution "multiPhATE") directory, and clone multiPhATE from the command line: 

```
$ git init

$ git clone https://github.com/carolzhou/multiPhATE2
```

(Complete instructions for using git and github can be found at http://help.github.com.)

Now, be sure that multiPhate.py and phate_runPipeline.py and associated files and directories are in your main execution "multiPhATE" directory. Check that the two subdirectories: PipelineInput/ and PipelineOutput/ are present (should already exist in the downloaded distribution). Place your phage genome fasta files (genome1.fasta, genome2.fasta, etc.) into the PipelineInput/ subdirectory. Place your configuration file (ie, your copy of sample.multiPhate.config) in the main execution directory (same level as multiPhate.py). A word of caution here:  it is always best to name your files and fasta contigs as strings lacking any spaces or special characters, as third-party codes over which we have no control may balk when encountering odd characters or spaces. I have attempted to make the multiPhATE code robust with respect to odd characters in fasta headers, but there is no guarantee. If you run into problems that might stem from this issue, you may run your genome fasta files through script cleanHeaders.py (find it in the Utility folder). 

You will need to acquire one or more of the databases listed below under SUPPORING DATABASES (Phantome and pVOGs are included in the multiPhATE distribution, so it is possible to begin with just those), and the 3rd party codes listed under SUPPORTING 3rd PARTY CODES. You will need to acquire at least one of the supported gene finders, but it is recommended to run as many of the four gene finders as is feasible so that the results can be more meaningfully compared. You will need to specifiy the locations of the supporting data sets and codes in the multiPhATE config file (see multiPhate.config), and you will need to locate your genome file(s) to the PipelineInput/ subdirectory. Once you have acquired the third-party codes and databases, you will be ready to configure the multiPhate.config file.

#### HOW TO WRITE A CONFIGURATION FILE

Summary:
Availability and locations of supporting databases and codes are to be specified in a configuration file. A sample configuration file is provided, called "sample.multiPhate.config". Make a copy of this file and rename it accordingly (eg., myGenomeSet_multiPhate.config). Hereafter we refer to this file as, multiPhate.config. The multiPhate.config file is configured according to established default parameters (just about everything turned off initially). Any of the parameters may be modified (switches turned on or off) by assigning 'true' or 'false'. It is suggested that you turn swithes off, then install each supporting gene finder and database in turn and test the pipeline.

Procedure:

1) At the command line, make a copy of the file, sample.multiPhate.config, and name it appropriately (hereafter referred to as 'multiPhate.config'):  `$ cp sample.multiPhate.config multiPhate.config`.  Then, edit your config file as described below.

2) List of Genomes:
For each genome to be processed, provide six lines under "Genome List:" and before "END of list":  for each genome, you need to list the genome number, the name of the genome fasta file, the genome type (typically 'phage', but could be 'bacteria'), the species, if known (no spaces), the name of the genome, and a name for the output directory to hold this genome's output files (again, no spaces), in that order. You can simply copy/paste the six lines provided as many times as needed, and fill in the information appropriate for each genome.

3) Processing Information:
You may configure the pipeline to perform gene finding only, or gene finding plus functional annotation. For example, you may want to examine the results of multiple gene finders before going forward with functional annotation. In order to configure phate to run gene finding only, set translate_only to 'true'; in this way, only gene-calling and translation (to peptide sequence) will be performed. If you set translate_only to 'false', then the pipeline will not stop at the translation step, but will proceed with functional annotation of the predicted genes (ie, blast and/or hmm). Normally the genetic_code should be set to '11', for prokaryotic.

4) Gene Callers:
The gene_caller option specifies which gene caller's results (ie, gene calls) will be used for subsequent functional annotation. The choices are:  'phanotate', 'genemarks', 'prodigal', or 'glimmer'.  To run a gene caller, you must have acquired that third-party code and installed it locally for use with multiPhATE. For each gene caller you wish to have run, set that caller's parameter to 'true'. In the usual case, you will want to specify gene_caller='phanotate' for annotation of phage genomes. You may also provide your own gene calls (referred to as "custom"), but your input file must conform to a modified (simplified) gff format that is recognized by PhATE. A custom gene-calls file must be named according to your genome's name (e.g., myGenome.custom) and placed in the PipelineInput/ directory. If you select more than one gene caller (or select one and provide one of your own), then PhATE will compute the following gene-call sets: superset (non-redundant list of all gene-caller gene calls), consensus (set of calls that are in agreement among at least two gene callers), and commoncore (calls that are in agreement among all callers). You may select as the gene-call set to be forwarded for annotation any of the gene callers or 'superset', 'consensus', or 'commoncore'. Only one gene-call set will be annotated in a given multiPhATE run.

5) Annotation:
Set to 'true' each blast or hmm process that you want to be run. Note that you must have acquired the associated database, and in the next section (Databases) you must configure the location of each database. You may also set the desired blast parameters. The blast_identity sets the minimum identity that will be considered; any blast result below that threshold will be ignored. The hit_count parameters will determine how many top hits will be reported. (Note that recent releases of blast+ suggest considering at least 5 hits, and will generate a warning message if you set the hit count to less than 5.) You may select hmm search using phmmer, jackhmmer, or hmmscan. Currently hmm searches are performed only using the protein fasta databases or the pVOG profile database (future releases of multiPhate are expected to support additional databases).

6) Databases:
For each database that you have in-house, specify the full path/filename. Note that you may need to prepare in advance all blast databases by running the "makeblastdb" utility (see instructions with blast+ code for how to do that). MultiPhate will only run with blast+; it does not support legacy blast. For instructions where to download the databases, see the SUPPORTING DATABASES section below. Note that KEGG is available by license. Note also that in some cases additional files are required. In this case, place the additional file(s) in the same directory as the associated blast database. For example, place the NCBI accession2taxid file in the same directory as your NCBI virus genome file (see below). If you are downloading datasets that you anticipate using specifically with multiPhATE, then it is suggested, for convenience, that you save them in the Databases/ folder in the multiPhATE distribution, but any database can be located anywhere on your local system; you need only indicate in the multiPhate.config file the full path/filename for each database. Remember, the pVOGs and Phantome data sets are included in the multiPhATE distribution in the Databases/ folder, but you will need to run makeblastdb to render the datasets blast-able (`$ makeblastdb -help`). If you will be running hmmscan, then you will need to format the pVOG database accordingly: 

```
$ cat VOG\*.hmm > pVOGsHmmProfilesDB.hmm

$ mv pVOGsHMMprofilesDB.hmm ../.

$ cd ..

$ hmmpress pVOGsHmmProfilesDB.hmm
```

7) Verbosity:
You may up- or down-regulate verbosity in the multiPhate.config file, under "# VERBOSITY". This includes an option to clean the (voluminous) raw blast and hmm search data from the output directories. It is suggested that clean_raw_data, and phate_progress be set to 'true'. The phate_warnings and phate_messages, when set to 'true', will generate voluminous output; set these to 'true' only when trouble-shooting the pipeline. 

8) See INSTALLATION AND SET-UP CHECKLIST below.


#### PIPELINE EXECUTION

Run the PhATE pipeline at the command line by passing your multiPhate.config file as an argument to the multiPhate.py pipeline driver script, as follows: `$ python multiPhate.py multiPhate.config`


#### SUPPORTING DATABASES

It is recommended that the user acquire as many of the following sequence databases and associated codes as is feasible, although none are actually required to run the code. (You may specify "translate_only='true'" to do gene finding then translation, and then stop at that point.) Databases are listed with at least one way of acquiring them, but there may be additional sources, and it is possible to substitute subsets of blast databases (e.g., a subset of the NCBI gene database in place of Refseq Gene). Running the script, getDBs.py (under Utility/), may assist with downloading and running makeblastdb for the databases. (Disclaimer: NCBI occasionally changes what what exactly they are offering on their ftp site, and getDBs.py might not have kept up with these changes.) Remember that NR is a very large database; nobody should have two copies of this database on disk. Script getDBs.py will allow you to select the specific databases you wish to download.

In some cases, you may prefer to download a database through a web interface, and in other cases you may use blast+ to download a database at the command line, if the database of interest is provided by NCBI. The latter method may be more convenient for larger data sets (eg, NCBI Refseq Protein). Blast+ includes a script (/bin/update_blastdb.pl), which can be invoked at the command line to automatically download a specified database. In order for blast+ to search a database, the database must first be converted to a blast-able object using the blast+ program, makeblastdb. Once you have installed blast+, you can query for help in using the program. For example, type at the command line: `$ makeblastdb -help`.

NCBI virus genomes - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/. If using a web browser, go to ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/. Then, select "Allow" to "Do you want to allow this page to open 'Finder'?". Connect as "Guest". Select the files to be downloaded:  viral.1.1.genomic.fna.gz, also 2.1, and 3.1. These can be dragged to your desktop, which should prompt the download.

NCBI-associated file:  accession2taxid file - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz. If using a web browser, go to ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. Then, select "Allow" to "Do you want to allow this page to open 'Finder'?". Connect as "Guest". Select the volumes to mount: "OK". This should download the zip file.

NCBI Refseq Protein - download using blast+: /bin/update_blastdb.pl refseq_protein

NCBI Refseq Gene - download using blast+: /bin/update_blastdb.pl refseqgene. This database contains primarily human sequences. To acquire bacterial gene sequences you may need to download them via http://ncbi.nlm.nih.gov/gene and process the data to generate a fasta data set. Support for doing this is not provided with the multiPhATE distribution.

NCBI Swissprot - download using blast+: /bin/update_blastdb.pl swissprot

NR - ftp://ftp.ncbi.nlm.nih.gov/nr/

KEGG virus subset - (available by license) http://www.kegg.jp/kegg/download/

KEGG associated files - T40000.pep, T40000.nuc, vg_enzyme.list, vg_genome.list, vg_ko.list, vg_ncbi-geneid.list, vg_ncbi-proteinid.list, vg_pfam.list, vg_rs.list, vg_tax.list, vg_uniprot.list

Phantome protein fasta sequences - http://www.phantome.org/Downloads/phage_proteins_nnnnnnnnn.fasta. (A version of Phantome is included in the multiPhATE distribution.)

pVOGs prepared database (pVOGs.faa) - included in PhATE distribution. This data set was derived by C. Zhou from the pVOGs fasta database. For use in PhATE, the sequence fasta headers have been modified to include the pVOG identifiers (all groups to which each sequence belongs). This re-formatting facilitates pVOG group identification and construction of the alignment-ready fasta files. Codes for reconstructing this modified data set are included in the PhATE distribution. Note that the pVOGs are not mutually exclusive, meaning that a sequence may have membership in more than one VOG group. The codes included in the phate distribution will combine identifiers that belong to a given sequence and list all the VOG identifiers in the fasta header. In this way, the pVOG fasta database used in PhATE is also non-redundant. See documentation in DatabasePrep/dbPrep_createPvogFastaFile.py for instructions how to update your local pVOGs data set for use in PhATE, but you can start with the pVOGs.faa file included in the PhATE distribution. Combine the pVOG fasta sequences into a single file and format for hmmscan profile search as follows:

```
$ cat VOG\*.hmm > pVOGsHmmProfilesDB.hmm
$ mv pVOGsHMMprofilesDB.hmm ../.
$ cd ..
$ hmmpress pVOGsHmmProfilesDB.hmm
```

For simplicity in configuring the locations of dependent databases in the multiPhate.config file, it is suggested that the above databases be placed in a directory structure as follows: 

```
Databases/
	KEGG/ 
	NCBI/ 
		Virus_Genome/ 
		Virus_Protein/ 
	NR/ 
	Phantome/ 
	Refseq/ 
		Protein/ 
		Gene/ 
	Swissprot/ 
	pVOGs/
```

You must specify in your multiPhate.config file the locations of the data sets that you will be using. Although it is recommended that you place your databases in the above directory structure, they can reside anywhere locally on disk, but in any case you must specify the full directory path/filename to a given resource in your multiPhate.config file.


#### SUPPORTING 3rd PARTY CODES

Note that some third-party codes are required for multiPhATE, but others are optional, as indicated below. Some of these codes can be installed in a Conda environment. Codes that can be installed via Conda are so indicated below. If using Conda, follow the instructions that occur at the bottom of this section. Otherwise, install these codes globally, following the instructions provided with each package from the source.

BioPython - https://biopython.org/wiki/Download (required) (conda)

EMBOSS package - https://sourceforge.net/directory/os:mac/?q=EMBOSS (required) (conda)

Blast+ https://ncbi.nlm.nih.gov/blast (optional) (conda)

GeneMarkS - http://exon.gatech.edu/Genemark/index.html (optional; available by license)

Glimmer - https://ccb.jhu.edu/software/glimmer/ Use Glimmer version 3. (optional) (conda)

Prodigal - https://github.com/hyattpd/Prodigal (optional) (conda)

PHANOTATE - A Python 3-compatible version of PHANOTATE is included in the multiPhATE distribution, in the ExternalCodes/ folder. Unzip PHANOTATE.zip, and follow instructions in the README. Future updates to PHANOTATE will be available at https://github.com/deprekate/PHANOTATE. (optional)

jackhmmer, phmmer, hmmscan - https://www.eddylab.org/software.html or http://www.hmmer.org/download.html. Download HMMER; hmm search codes are included in this package (optional) (conda - hmmer) 

tRNAscan-SE - https://www.eddylab.org/software.html - select tRNAscan-SE download link (conda)

Third-party codes should be installed globally whenever possible. However, it is recommended that PHANOTATE be installed under the ExternalCodes/ subdirectory in the execution/working directory. (The ExternalCodes/ subdirectory should already exist in the multiPhATE distribution, with PHANOTATE for Python 3.x in that location.)

#### CONDA INSTALLATION

If you prefer to run multiPhATE in a Conda environment, here are some tips for how to set it up. 

1) First, download and install miniconda3 for Python 3.7 (https://conda.io/en/latest/miniconda.html). For more information about Conda, see https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html.  

2) Check that conda is working:  `$ conda --version`
 
If conda is not recognized, then you may need to  switch to bash shell: `$ bash`  
(and then try again)

3) Add the following channels:  

```bash
$ conda config --add channels defaults 
$ conda config --add channels conda-forge 
$ conda config --add channels bioconda 
```

Note: bioconda is supported on Linux and Mac operating systems, but so far not on Windows.

4) Create an environment for using multiPhATE; let's call it "multiphate":  `$ conda create --name multiphate`

5) Activate that environment:  `$ source activate multiphate`

6) Install conda packages within that environment:  `$ conda install python=3`

Repeat for each of biopython, emboss, blast, glimmer, prodigal, hmmer, trnascan-se.  

7) When running multiPhATE within your multiphate Conda environment, the pipeline will use the version of python and the third party codes installed within the multiphate environment, so there should be no clashes with other versions of these packages that may be installed elsewhere on your system. When you are finished running multiPhATE, you may exit from the multiphate Conda environment:  $ source deactivate

Note that genemarks and phanotate are not available as conda packages, so these programs, as well as the dependent databases, all need to be acquired/installed manually in any case.

#### PHATE PIPELINE OUTPUT FILES 

In the user-specified output directory (eg, myGenomeDir/), under PipelineOutput/, the following files will be written:

gene-call outputs from each of the gene callers that was run, including a gff-formatted output file

gene.fnt and protein.faa fasta files generated using the designated preferred gene finder

CGC_results.txt - side-by-side comparison of all gene finder results (if at least two were run) plus comparison statistics

cgc.gff - a superset of gene calls, each tagged with the gene caller(s) that made the call

phate_sequenceAnnotation_main.out - tabbed integrated annotation results (also written in gff format)

BLAST/ directory - raw blast results for genome (under Genome/) and proteins (under Protein/) (if CLEAN_RAW_DATA switch in multiPhate.config file is 'false')

pVOG groupings comprising alignment-ready fasta files (under BLAST/Protein/) based on best blast hits

HMM/ directory - raw hmm search results for protein sequences (under Protein/) (if CLEAN_RAW_DATA switch in phate_runPipeline.py is 'False') This includes the hmm-formatted outputs (tbl and dom) and the raw alignments (std), which are voluminous.

pVOG grouping comprising alignment-ready fasta files (under HMM/Protein/) based on best hmm hits

log files capturing details of the processing and time stamps

The auto-generated myGenomeName_phate.config file, to record exactly how you configured the pipeline for the current run (genome).

#### INSTALLATION AND SET-UP CHECKLIST

* Have you installed multiPhATE either by downloading the zip file from https://github.com/carolzhou/multiPhATE.git or cloning the repository?
* Have you installed the databases you would like to use with multiPhATE? Recall that Phantome and pVOGs are included under Databases/, and that script getDBs.py can assist you in downloading and preparing the databases.
* Have you run makeblastdb on each of your local databases (if this has not already been done by getDBs.py)?
* Have you acquired the dependent codes, either by installing from the provided web addresses or by installing within a Conda environment?
* Have you created a copy of sample.multiPhate.config and configured it?
	- added meta-data for each genome to be processed
	- configured "translate_only" to do genecalling only ('true') or to do genecalling followed by blast/hmm analyses ('false')
	- selected gene callers to be run, and specified the gene caller to use for annotation (preferred caller)?
	- specified the blast/hmm analyses to be performed
	- specified the locations/names of the databases you have locally on your system
	- modified the verbosity (optional)
* We recommend stepwise testing to be sure all components have been correctly installed and specified.
* Feel free to post issues and suggestions regarding multiPhATE on our github project page: https://github.com/carolzhou/multiPhATE2.git. Select the 'Issues' tab.

#### RUNNING PHATE AS AN "EMBARASSINGLY PARALLEL" CODE

There are three ways to perform parallel processing with multiPhATE. 1) The code can be run using pthread by specifying the number of threads in the configuration file (threads=''), or by specifying 'ALL', to use all available pthreads on your system. 2) Multiple instances of multiPhATE can be distributed across cores of a high-performance computing machine by specifying HPC='true' in the configuration file. Note that the user will need to write scripts specific to the hardware on which multiPhATE is to be run. 3) Blast+ allows the user to specify the number of threads for running blast processing. Specify the number of blast threads in the configuration file: blast_threads=''. 

The MultiPhATE2 Pipeline avoids clashes in writing results; outputs for each genome are written to user-specified output subdirectories (specified in your multiPhate.config file).

#### FURTHER RECOMMENDATIONS

PhATE was originally developed for anotating phage genome sequences. However, PhATE may also be useful in helping to identify phage genes within bacterial genomes (i.e., prophage). Thus, the user has the option of running multiple gene callers for bacterial genome sequence (GeneMarkS, Glimmer, Prodigal) and a new gene finder specifically for phage (PHANOTATE). The calls from any two or more of these callers are compared by the PhATE/CGC code so that the user can examine the calls that agree or disagree among the callers, and then run PhATE again, selecting the caller of choice for annotating the sequence.

Although most supporting databases for multiPhATE are phage- or virus-centric, NR and Refseq Protein are included in order to help identify genes/functions that are not known in the virus/phage gene data sets. However, PhATE is not intended to be the sole source for anntotation of bacterial genome sequences, as PhATE is tailored for identification of genes and functions in phage.

The NR database has grown enormously large. It is recommended to use a smaller database, such as Refseq Protein instead of NR. Furthermore, annotating with NR will add greatly to the time required for processing a genome through PhATE. Therefore, it is recommended that NR be turned off ('false') until one desires to preform a full/final annotation of the genome of interest, if using NR.

Because the behavior of 3rd party codes can sometimes be unpredictable, it is recommended that the user replace spaces and unusual characters in their fasta headers with the underscore character.

#### CAUTIONS 

1) We have detected an issue with PHANOTATE/fastpath. There appears to be a system-dependent problem in compiling fastpath, which disables PHANOTATE. This error occurs on some systems but not others. If running multiPhate.py yields an error message "IndexError: list index out of range" in phate_geneCallPhage.py, line 216, please disable phanotate in your multiPhate.config file (phanotate_calls='false', and gene_caller='\<your choice other than phanotate\>'), and continue evaluating multiPhATE without PHANOTATE for the time being. When this issue has been resolved we will upload a patched code and will delete this message from the README. Thank you.

2) multiPhATE2 has not been tested on Windows. 

#### PLANNED FURTHER DEVELOPMENT

Plans include adding additional functionality to the pipeline:
1) Additional hmm databases.
2) Computing data files for a reference genome comprising each gene's closest homologs in each of the other specified genomes.
3) Creating hmm profiles (via hmmbuild) for each of a reference genome's genes + closest homologs.
4) Implementing optional "custom" databases for blasting or doing hmm analysis. (Actually, you can do custom blast now by substituting your custom nucleotide blast database for the refseq_gene_blast process, or substituting your custom protein blast database for one of the protein processes (swissprot, nr, refseq_protein, or ncbi_virus_protein), by specifying the location of your custom database in the multiPhate.config file and turning the corresponding annotation process to 'true'.)

Feel free to report bugs or problems, or to suggest future improvements, by posting an issue on the project github page (click on the Issues tab), or by emailing the developers at:  zhou4@llnl.gov. Thank you for using multiPhATE.

#### PUBLICATION

If you use multiPhATE in your research, kindly reference our paper:  "multiPhATE: bioinformatics pipeline for functional annotation of phage isolates", by Carol E Zhou, Stephanie A Malfatti, Jeffrey A Kimbrel, Casandra W Philipson, Katelyn E McNair, Theron C Hamilton, Robert A Edwards, and Brian E Souza; Bioinformatics, 2019, https://doi.org/10.1093/bioinformatics/btz258. If you run PHANOTATE "under the hood", please also reference, "PHANOTATE: a novel approach to gene identification in phage genomes", by Katelyn McNair, Carol Zhou, Elizabeth A Dinsdale, Brian Souza, and Robert A Edwards, Bioinformatics, 2019, https://doi.org/10.1093/bioinformatics/btz265. 

multiPhATE2 v.1.5
