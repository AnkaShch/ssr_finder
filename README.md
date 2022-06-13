# SSR Finder
## Introduction
ssr_finder is a tool for searching for short tandem repeats for a given motive in FASTA sequence files.

## Installation
ssr_finder can be installed from the source code:
```
# Download the git repo
$ git clone https://github.com/AnkaShch/ssr_finder.git
```
## Usage
The help message and available options can be accessed using
```
$ python ./__main__.py --help
$ python ./--main__.py -h
```
which gives the following output
```
usage: __main__.py [-h] [--motif MOTIF] [--motif_file MOTIF_FILE]
                   [--output_prefix OUTPUT_PREFIX] [--region REGION]
                   [--cores CORES] [--valid_regex] [--distance DISTANCE]
                   [--number NUMBER] [--search_strand {+1,-1,both}]
                   [--no_motif_overlaps]
                   [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   fastaFile

positional arguments:
  fastaFile             Path to FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  --motif MOTIF, -m MOTIF
                        A degenerate sequence motif. Can be specified multiple
                        times.
  --motif_file MOTIF_FILE
                        A file containing motifs to search file. Can be
                        specified multiple times.
  --output_prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Prefix of generated bed files. Default is
                        detected_ssrs
  --region REGION, -r REGION
                        Defines region to search for sites. Use 'contig:start-
                        end' for regions, or 'contig' for whole contig. If no
                        regions are specified, the entire FASTA file will be
                        searched! Starts expected to be *0 start* start and
                        Ends *1 start*
  --cores CORES, -p CORES
                        Run search on multiple contigs / strands
                        simultaneously
  --valid_regex         Query is valid regex. *WILL NOT* reverse complement.
                        Specify sequence and strand with careful
                        consideration.
  --distance DISTANCE, -d DISTANCE
                        The maximum allowable distance between motifs in bp,
                        within which it is assumed that it all falls into the
                        repeat region. Default is 0.
  --number NUMBER, -n NUMBER
                        The minimum number of motives in one SSR. Default is
                        2.
  --search_strand {+1,-1,both}, -s {+1,-1,both}
                        Default searches both strands, but can be set to only
                        one.
  --no_motif_overlaps   If this option is set it turns off overlapping motif
                        matches.
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}

__main__.py v0.5
```
The details of each option are given below:

**`fastaFile`**

**Expects:** _FASTA file_

**Default:** _None_

It is position argument with path to FASTA file. It has not specific flag.

**`--output_prefix or -o`**

**Expects:** _STRING (to be used as prefix of filename)_

**Default:** _detected_ssrs.bed and detected_ssrs_full.bed_

If this option is not provided, default outputs filename will be `detected_ssrs.bed` and `detected_ssrs_full.bed`. 
Please note that even in the case of no identified SSRs, output files is still created (therefore overwriting any previous file of the same name) but with no content in the file.
The Tool creates 2 files in bed format at the output.
One of its with standard data (more details below), the second (with the suffix `_full.bed`) contains additional information about the found repeats.

\\\ TBD
