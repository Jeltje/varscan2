# Varscan2

**This repository contains code to create a docker implementation of the Varscan2.4.2 copynumber variation (CNV) caller.**

Varscan2 was developed by Dan Koboldt (see References below). It can be used to detect copy number variation (CNV) in sample pairs, usually exomes from a tumor and control from one patient.

The Varscan2 executable (https://github.com/dkoboldt/varscan.git) combines several tools. It is meant to be run in a pipeline, during which different tools are called in sequence. For details on Varscan, see http://dkoboldt.github.io/varscan/


This repository ONLY contains a pipeline for Varscan2 **copynumber variation**. If you want to run other Varscan tools, please use Varscan2 directly. This docker container contains a wrapper script that uses Varscan tools and other programs *with specific parameters*. These may not be the perfect parameters for your particular samples. See below for the full list of pipeline steps.

**Inputs** to the program are a tumor/control pair of BAM files and several [bed format](https://genome.ucsc.edu/FAQ/FAQformat#format1) helper files (see below).
**Output** is a file with chromosome segments that are scored for amplification or deletion.

To get per-gene output, these scores must be mapped to an annotation, for example using [this program] (https://github.com/Jeltje/cnvtogenes)

## The code

The Varscan wrapper script runs the following:

1. samtools flagstat on each bam file
2. samtools mpileup on both bam files
3. Determine unique mapped read ratio
4. Varscan copynumber
5. Remove low coverage regions
6. Varscan copyCaller
7. Calculate median for recentering
8. Varscan copyCaller recenter
9. Separate chromosome arms
10. DNAcopy (CBS)
11. Merge chromosome arms

The chromosome arms are separated before the Circular Binary Segmentation (CBS) step to avoid making calls across centromeres.

## Getting the docker container

The latest Varscan docker image can be downloaded directly from quay.io using
`docker pull quay.io/jeltje/varscan`

Alternatively, you can build from the github repo:
```
git clone https://github.com/jeltje/varscan2.git
cd varscan2
docker build -t jeltje/varscan .
```

## Running the docker container

For details on running docker containers in general, see the excellent tutorial at https://docs.docker.com/

To see a usage statement, run

```
docker run jeltje/varscan -h
```

### Example input:

```
docker run  -v /path/to/input/files:/data jeltje/varscan -c normal.bam -t  tumor.bam -q sampleid -i genome.fa -b centromeres.bed -w targets.bed -s tmpdir > varscan.cnv

```

where

`normal.bam` and `tumor.bam`    are BAM format files of exome reads aligned to the genome. 

`sampleid` is an identifier for the patient. This will be used in the output.

`genome.fa` is a fasta file containing the genome that was used to create the BAM files. A samtools indexed `.fai` file must be present in the same directory as this file (for details see Other Considerations, below)

`centromeres.bed` is a [bed format file](https://genome.ucsc.edu/FAQ/FAQformat#format1) containing centromere locations. This list is used to remove centromeres from the CBS calls.

`targets.bed` is a list of exome targets in bed format. This is used as a 'whitelist' of genome regions so that off target alignments will not be used for analysis

`tmpdir` is a directory for temporary output files. If you set option -d, these files will be kept

Keep in mind that all these file locations must be with respect to your `/path/to/input/files`.

Centromeres for hg19 are provided ind the `/data` directory

>       You can find centromere locations for genomes via
>       http://genome.ucsc.edu/cgi-bin/hgTables
>       Using the following selections:
>       - group: Mapping and Sequencing
>       - track:gap
>       -       filter - goes to new page, look for 'type does match' and type centromere, submit
>       -       output format: bed
>       Submit, on the next page just press Get Bed


## Output

Output is written to `STDOUT` and uses the following format:
```
sampleID    chrom    loc.start    loc.end    num.mark    seg.mean

```

To get amplified or deleted segments from this file, a threshold must be applied. This is often set to `0.25/-0.25`,
and with a minimum number of 10 markers per segment.


## Other considerations

Tumor and control really must be from the same patient and processed in the same experiment. Batch effects are strong in exome experiments and using the wrong control renders Varscan output meaningless.

To index a genome, run
```
	samtools faidx <genome.fa>
```
This creates an index named genome.fa.fai
The genome and the index must be in the same directory, and the genome file (not the index) is the input to run_varscan

The whitelist is a bed format file with the exome targets used in the experiment. It ensures that Varscan only uses target regions for its analysis and not any off target read matches. It is important to use the real list of exome targets. For meaningful results do not use a generic list.


## References

Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. 
VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. 
Genome Res. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111.
