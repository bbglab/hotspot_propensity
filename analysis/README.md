# Analysis

This directory contains different analyses carried out in the manuscript. 

Note that part of the analysis has been run in a high computing cluster. For such cases, Python scripts are provided together with Jupyter Notebooks that prepare job submission. 

## Input data

The input data to run these scripts is referenced in the manuscript. You will need to download the data from the original sources to run different parts of the analysis. Take also into account that some parts of the analysis are linked to one another, therefore have to be run in a given order as shown below. 

## Content

#### 1) Cancer driver genes and regions

Annotations of cancer driver genes and their associated coding and non-coding genomic regions. 

Folder: cancerdrivers

#### 2) Mappable genome

Code to generate mappable genome annotations and compute their trinucleotide composition.

Folder: mappable\_genome

#### 3) Pre-processing somatic mutations from cancer cohorts

This directory includes:
- code to filter samples and somatic mutations from cancer cohorts. 
- code to merge cohorts into cancer types
- analysis of shared mutations between cancer types
- sample annotations

Folder: processing\_cancer\_cohorts

#### 4) Hotspots across cancer genomes

Code to run HotspotFinder across cancer types.

Folder: hotspots\_cancer\_types

#### 5) Genomic bins

Code to generate 1 Mbp and 500-10 Kbp bins of the mappable genome. Includes code to generate regional genomic annotations of trinucleotide genome composition, mutation rates, hotspot rates, chromatin accessibility, replication timming, and gene expression. 

Folder: genomic\_bins

#### 6) Extraction of mutational signatures from cancer genomes

This folder contains the code for the generation of mutational profiles with SigProfilerMatrixGenerator, extraction mutational signatures with SigProfilerExtractor, assignment of mutations to signatures, and assignment of hotspots to signatures. 

Folder: mutational\_signatures

#### 7) Hotspot propensity 

This folder contains the code to calculate the observed hotspot propensity across 14 signatures in 7 cancer types. 

Folder: hotspot\_propensity

#### 8) Expected hotspot propensity 

This folder contains the code to generate the expected estimates of hotspot propensity across cancer genomes, including 1 Mbp and 500-10 Kbp bins models for different cancer type and signature pairs, as well as methylation aware models for SBS1. 

Folder: theoretical\_models

#### 9) CTCF

Analysis of SBS17 hotspots and tissue-matched CTCF binding sites. 

Folder: ctcf

#### 10) Methylation

Analysis of SBS1 hotspots and tissue-matched CpG methylation. 

Folder: methylation

