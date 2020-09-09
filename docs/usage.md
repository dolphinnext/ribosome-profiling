# dolphinnext/ribosome-profiling: Usage

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run dolphinnext/ribosome-profiling -profile docker --DOWNDIR /path/to/save/ribosome-profiling --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build human_hg38_gencode_v30
```

If you're running for the first time, you need to enable `--run_checkAndBuild` paramater as follows:

```bash
nextflow run dolphinnext/ribosome-profiling -profile docker --DOWNDIR /path/to/save/ribosome-profiling --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build human_hg38_gencode_v30 --run_checkAndBuild 'yes'
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull dolpinnext/ribosome-profiling
```

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from Dockerhub: [`dolphinnext/riboseq`](http://hub.docker.com/r/dolphinnext/riboseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character


### `--mate`
You can specify the single-end data by entering mate parameter as 'single'. For example:

```bash
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```


## Reference genomes

### `--genome_build` 
To run the pipeline, you must specify which to use with the `--genome_build` flag.

List of genomes that are supported are:

* Human
  * `--genome_build human_hg38_gencode_v30`
  * `--genome_build mouse_mm10_gencode_m25`

Note: For new genome requests, please send e-mail to UMMS-Biocore(biocore@umassmed.edu).


## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human hg19 reference genome.

You can choose multiple aligner to compare their results by enabling/disabling following parameters:
```bash
To enable STAR    : `--run_STAR yes`
To disable STAR   : `--run_STAR no`
```

## Adapter Removal
If specific Adapter Removal is required, you can enable trimmomatic and enter the adapter sequence. 

```bash
To enable adapter_removal: 
--run_Adapter_Removal "yes"

--Adapter_Trimmer_Quality_Module_Adapter_Removal.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence [string]
# You can enter a single sequence or multiple sequences in different lines. Reverse sequences will not be removed.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length [int @default:10]
# Specifies the minimum length of reads to be kept

--Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches [int @default:1]
# Specifies the maximum mismatch count which will still allow a full match to be performed

--Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold [int @default:30]
# Specifies how accurate the match between the two -adapter ligated- reads must be for PE palindrome read alignment

--Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold [int @default:5]
# Specifies how accurate the match between any adapter etc. sequence must be against a read.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped [@options:"yes","no" @default:"yes"]
# Discard_non_clipped sequences (keep only sequences which contained the adapter)
```

## Trimmer
Optianally, you can trim your reads by defining trimming lenghts as shown at below: 

```bash

--run_Trimmer [@options:"yes","no" @default:"no"]
# Enables Trimmer by setting this parameter as "yes"

--Adapter_Trimmer_Quality_Module_Trimmer.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

For Single End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "single"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime [int]

For Paired End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "pair"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2 [int]
```

## Quality Filtering
Optianally, you can trim your reads based on their quality. Trimmomatic works on both paired-end and single ended data. Alternatively fastx option (fastx_toolkit) could be used for single reads. 

```bash
To use Trimmomatic  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "trimmomatic"

--Adapter_Trimmer_Quality_Module_Quality_Filtering.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size [int @default:10]
# Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold (=required_quality).

--Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming [int @default:15]
# Specifies the average quality required for window trimming approach

--Adapter_Trimmer_Quality_Module_Quality_Filtering.leading [int @default:5]
# Cut bases off the start of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing [int @default:5]
# Cut bases off the end of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen [int @default:36]
# Specifies the minimum length of reads to be kept
```

```bash
To use fastx_toolkit  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "fastx"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality [int @default:20]
# Minimum quality score to keep reads

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent [int @default:100]
# Minimum percent of bases that must have entered minQuality
```

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). 

Note - you can use this to override pipeline defaults.