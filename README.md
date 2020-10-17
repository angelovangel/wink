# WINK

What's In my Nanopore reads, with Kraken2, in real-time

## Description

WINK is a platform for real-time phylogenetic classification and species quantification for Nanopore sequencing data, based on [kraken2]() and [bracken](https://ccb.jhu.edu/software/bracken/). It can be used both in real-time (monitor a specified folder for new reads, e.g. `fastq_pass` and continuously update results) and post-run (collect all the reads and perform analysis).<sup>[1](#footnote1)</sup> The software consists of two parts - a nextflow pipeline (can be executed on its own) and a graphical user interface (a Shiny app) which collects the output of the nextflow pipeline and diplays it as a dashboard page.


## Running the pipeline and explanation of the results

WINK can be run either via the Shiny app (in a browser, no command line needed), or by executing the nextflow pipeline from the command line.

The input for the pipeline is one of:

- the output folder of the MinKNOW software (usually this is the `fastq_pass` folder) or 
- any other folder where the basecalled data (fastq files) accumulates during a run.
*The run may be barcoded or not*

The results can be monitored in real time via the app. In addition, the nextflow pipeline outputs in real-time (all in `results-wink`):

- `latest-fastq` - where the `latest` fastq files are collected, one file per sample. By default, the MinKNOW software writes 4k reads per fastq file and typically one sample can generate many such fastq files. The WINK pipeline merges all the available fastq files belonging to one sample, i.e.

```bash
data
├── barcode01
│   ├── PAE58908_pass_barcode01_d2c5c063_0.fastq
│   ├── PAE58908_pass_barcode01_d2c5c063_1.fastq
│   ├── PAE58908_pass_barcode01_d2c5c063_2.fastq
```

becomes `barcode01.fastq`. The `latest-fastq` is continously updated as new reads are generated.

- `latest-stats` - tables (one per sample/barcode) with the following columns:
`format type num_seqs sum_len min_len avg_len max_len Q1 Q2 Q3 sum_gap N50 Q20(%) Q30(%)`. These are also continuously updated.

- `latest-bracken`- a table with the bracken quantification of species abundance in the samples with the following columns:
`file name taxonomy_id kraken_assigned_reads new_est_reads freq`

## Under the hood

The pipeline is built with [nextflow](https://www.nextflow.io/) and [Shiny](), using some built-in nextflow functions to watch for new reads in the input folder. As new reads are generated, their phylogenetic assignment is performed with kraken2 and the relative species composition is determined with bracken. In parallel, various statistics about the fastq reads are collected and updated during the run. WINK runs in a docker container, so there are no external requirements.<sup>[2](#footnote2)</sup>


<a name="footnote1">1</a>: The reads generated *before* the pipeline is started are also included in the analysis by changing their timestamps. Take care if you rely on this information for other purposes.

<a name="footnote2">2</a>: Certainly will not run on Windows. MacOS and Linux are fine.