Ribosome Profiling pipeline is for processing Ribo-seq data. It accepts an adapter sequence for adapter removal and a fastq file(s) from Ribo-Seq experiment. 
It uses STAR aligner to ncRNA removal and rRNA mapping. Main features:
* To measure readthrough of stop codons genome-wide; average gene (metagene) analysis is performed by aligning all transcripts at their annotated stop codons and calculating normalized ribosome densities in this window.
* To evaluate Stop codon readthrough on a per transcript basis; Ribosome ReadThrough Score (RRTS) is calculated which is the density of ribosomes in the region 
of the 3′UTR between the Normal Termination Codons and the first in-frame 3′TC, and divided this value by the density of ribosomes in the CDS for every annotated 
transcript.

This pipeline adapted from following study: <a class="link-underline" href="https://elifesciences.org/articles/52611" target="_blank">Paper</a>  and <a class="link-underline" href="https://github.com/jrw24/G418_readthrough" target="_blank">Code</a> 


##### Citation:

* If you use DolphinNext in your research, please cite: 
Yukselen, O., Turkyilmaz, O., Ozturk, A.R. et al. DolphinNext: a distributed data processing platform for high throughput genomics. BMC Genomics 21, 310 (2020). https://doi.org/10.1186/s12864-020-6714-x
* Wangen, J.R., and Green, R. (2020). Stop codon context influences genome-wide stimulation of termination codon readthrough by aminoglycosides. Elife 9. 2020;9:e52611.

##### Steps:

1. For Quality Control, we use FastQC to create qc outputs. 
2. Adapter Removal: Adapter removal is performed by trimmomatic.
3. There are optional read quality filtering (trimmomatic) and read quality trimming (trimmomatic) processes available after adapter removal. 
5. ncRNA removal: noncoding sequences (Mt_rRNA, Mt_tRNA, rRNA, miRNA, scRNA, scaRNA, snoRNA, snRNA, sRNA, vaultRNA) removed by STAR.
6. STAR genome alignment: Remaining reads were mapped to the genome using STAR.
7. Separate density files for each read length around 15-40 nuc. created.
8. Metagene around the first inframe start/stop codon is created.
9. By using density files, normalized reads around the start/stop codon plotted.
10. Codon occupancies relative to the control sample plotted.
11. Read size distributions for each region of an mRNA plotted.
12. Transcripts sorted by stop codon identity and measured RRTS values for all transcripts.

##### Inputs:

* **Reads**: Specify the location of your input FastQ file. <a class="link-underline" href="https://dolphinnext.readthedocs.io/en/latest/dolphinNext/quick.html#adding-files" target="_blank">Need Help?</a>
* **Adapter Sequence**: Please enter the adapter sequence(s) in the settings of `run_Adaper_Removal`.
* Settings of `run_riboseq_workflow`: 
	- `sample_order` (optional): You can overwrite the default order of the samples in the figures by entering a new set of 'comma-separated' name of the samples. e.g. `control_rep1, control_rep2`
	- `amino_acid_list` (required): Please enter comma-separated list of amino acids that are going to be highlighted in figure 2S3B.
	- `color_code_list` (optional): 
	- `control_group_name`  (required): Control group name for figures  e.g. `control`
	- `control_group` (required): Comma-separated list of sample names e.g. `control_rep1, control_rep2`
	- `treatment_group_name`  (required): Treatment group name for figures (e.g. for first group: `treatment1` and for second group click add button and enter:`treatment2`)
	- `treatment_group` (required): Comma-separated list of samples (e.g. for first group: `treat1.rep1, treat1.rep2` and for second group enter: `treat2.rep1, treat2.rep2`)

##### Program Versions:
  - fastqc=0.11.8
  - star=2.6.1d 
  - samtools=1.6
  - multiqc=1.7
  - trimmomatic=0.39
  - bedtools=2.27.1
  - pandas=0.22.0
  - ucsc-fatotwobit=377
  - argparse=1.4.0
  - pysam=0.13
  - scipy=1.0.1
  - statsmodels=0.9.0
  - twobitreader=3.1.4
  - pathos=0.2.1
  - matplotlib=2.2.2
  - seaborn=0.9.0
  - logomaker=0.8
  - scikit-learn=0.20.3
  - r-ggplot2=2.2.1
  - r-plyr=1.8.4
  - r-reshape2=1.4.3
  - r-scales=0.5.0
  - xtail=1.1.5

##### Run through DolphinNext User Interface:

To start using the dolphinnext/ribosome-profiling pipeline please go to <a class="link-underline" href="https://dolphinnext.umassmed.edu/index.php?np=1&id=688" target="_blank">DolphinNext Web page</a> and click run button.

##### Run through Command Line:

To install and start using the dolphinnext/ribosome-profiling pipeline by using command line, please follow these steps: <a class="link-underline" href="https://github.com/dolphinnext/ribosome-profiling/blob/1.0/docs/local.md" target="_blank">Installation</a> .