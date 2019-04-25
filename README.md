# RNAseq-NF pipeline 

A basic pipeline for quantification of genomic features from short read data
implemented with Nextflow adapted for PILM Bioinformatics.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/pilm-bioinformatics/pipelines-nf-fastrnaseq.svg?branch=master)](https://travis-ci.org/pilm-bioinformatics/pipelines-nf-fastrnaseq)

## Requirements 

* Unix-like operating system (Linux, macOS, etc)
* Java 8 

## Quickstart 

Full tutorial lives at our [knowledgeBase](https://pilm-bioinformatics.github.io/knowledgebase/analysis/rnaseq/).

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).

2. Install Nextflow (version 0.24.x or higher):
      
        curl -s https://get.nextflow.io | bash

3. Assuming you have conda installed:

	conda -n create -y nf-fastrnaseq nextflow salmon multiqc fastqc
	conda activate nf-fastrnaseq

	
You can see an example report at the following [link](http://multiqc.info/examples/rna-seq/multiqc_report.html).	
	
Note: the very first time you execute it, it will take a few minutes to download the pipeline 
from this GitHub repository and the the associated Docker images needed to execute the pipeline.  


## Cluster support

RNASeq-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an 
abstraction between the pipeline functional logic and the underlying processing system.

This pipelines suppors `local` and `slurm` at eofe MIT cluster.

## Components 

RNASeq-NF uses the following software components and tools: 

* [Salmon](https://combine-lab.github.io/salmon/) 0.8.2
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Multiqc](https://multiqc.info) 1.0

