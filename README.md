# RHEA: Reference-free Heterogeneity and Evolution in Assembly graphs

### Description

Rhea is a software used to detect structural variants (SVs) between steps in long-read metagenomic series data. Current 
SVs detected include: insertions, deletions, mutations, and tandem duplications.

### Demo

Detect SVs between two samples
```bash
python rhea.py ..
```

Expected output:

Expected run time:

### Installation

##### Dependencies



### Parameters

| Command	| Default	| Description	|
| :-------  | :----- | :-------- | 
|input | (required)	| paths to input metagenome sequences or alignments|
|--output-dir, -o | rhea_results | path to output directory |
|--type | nano-raw | type of reads for MetaFlye graph construction ['pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw', 'nano-corr', 'nano-hq'] |
|--input-graph | (generated) | path to graph if alignments are provided|
|--bp-table | (generated)	| path to bp counts per sample if alignments are provided|
|--node-std | 1	| number of standard deviations away from median to call indels and mutations|
|--edge-lfc-thresh | 1	| number of lfc increase to call duplications |
|--raw-diff | FALSE | set to true if no normalization for bp count between samples is desired |
|--flye-exec | flye	| path to flye executable |
|--minigraph-exec | minigraph | path to minigraph executable |
|--threads, -t | 3| number of threads utilized by flye & minigraph|



### System Requirements

### Rhea Manuscript

