# SSR Finder
## Introduction
ssr_finder is a tool for searching for short tandem repeats for a given motive in FASTA sequence files. Our algorithm is based on the idea of the [motif_scaper](https://github.com/RobersonLab/motif_scraper) tool.

ssr_finder is crossplatformed tool. It is compatible with python 3 (tested on 3.8).

## Installation
ssr_finder can be installed from the source code:
```
# Download the git repo
$ git clone https://github.com/AnkaShch/ssr_finder.git
```
Install the necessary dependencies:
```
$ pip install pyfaidx==0.7.0
$ pip install regex>=2016.01.10
```

## Usage
Use a simple command to test the tool:

```
python __main__.py sample_data\data.fasta -m TTAGGG -o test_TTAGGG
```

The file `test_TTAGGG.bed` should contain the following:

```
#Contig	Start	End	Motif	Length	Strand	Number of motifs	Sequence
chr2	2	158	TTAGGG	156	+	26	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	163	181	TTAGGG	18	+	3	TTAGGGTTAGGGTTAGGG
chr2	182	206	TTAGGG	24	+	4	TTAGGGTTAGGGTTAGGGTTAGGG
chr2	209	263	TTAGGG	54	+	9	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	264	288	TTAGGG	24	+	4	TTAGGGTTAGGGTTAGGGTTAGGG
chr2	296	314	TTAGGG	18	+	3	TTAGGGTTAGGGTTAGGG
chr2	315	363	TTAGGG	48	+	8	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	364	418	TTAGGG	54	+	9	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	420	480	TTAGGG	60	+	10	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	485	545	TTAGGG	60	+	10	TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
chr2	554	566	TTAGGG	12	+	2	TTAGGGTTAGGG
```

The help message and available options can be accessed using
```
$ cd ssr_finder
$ python ./__main__.py --help
```
which gives the following output
```
usage: __main__.py [-h] [--motif MOTIF] [--motif_file MOTIF_FILE] [--output_prefix OUTPUT_PREFIX] [--region REGION] [--cores CORES] [--valid_regex] [--distance DISTANCE] [--number NUMBER] [--strand {+,-,both}] [--motif_overlaps]
                   [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   fastaFile

positional arguments:
  fastaFile             Path to FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  --motif MOTIF, -m MOTIF
                        A degenerate sequence motif. Can be specified multiple times.
  --motif_file MOTIF_FILE
                        A file containing motifs to search file. Can be specified multiple times.
  --output_prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Prefix of generated bed files. Default is detected_ssrs
  --region REGION, -r REGION
                        Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire FASTA file will be searched! 0-based format, start is included,
                        end in not included
  --cores CORES, -p CORES
                        Run search on multiple contigs / strands simultaneously
  --valid_regex         Query is valid regex. *WILL NOT* reverse complement. Specify sequence and strand with careful consideration.
  --distance DISTANCE, -d DISTANCE
                        The maximum allowable distance between motifs in bp, within which it is assumed that it all falls into the repeat region. Default is 0.
  --number NUMBER, -n NUMBER
                        The minimum number of motives in one SSR. Default is 2.
  --strand {+,-,both}, -s {+,-,both}
                        Default searches both strands, but can be set to only one.
  --motif_overlaps      If this option is set it turns on overlapping motif matches.
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}

__main__.py v1.0

```
The details of each option are given below:

**`fastaFile`**

**Expects:** _FASTA file_

**Default:** _None_

It is position argument with path to FASTA file. It has not specific flag.

**`--motif or -m`**

**Expects:** _STRING_

**Default:** _None_

This argument sets the motif for the search. Can be set multiple times.

The motif can be set in several ways:
1. the exact sequence of nucleotides: `-m TTAGGG`
2. using degenerate nucleotide codes:  `-m TTWGGG` for searching motifs TTTGGG and TTAGGG within one repeat
3. using python-style regular expression: `-m TTAG{2,3}` for searching motifs TTAGG and TTAGGG within one repeat

Important: regex string must be without spaces and this format requires flags `--valid_regex -s +1`.

**`--motif_file`**

**Expects:** _FILE_

**Default:** _None_

A file with a list of search motifs. Each motif is set on a separate line. Can be set multiple times.

**`--output_prefix or -o`**

**Expects:** _STRING (to be used as prefix of filename)_

**Default:** _detected_ssrs_

If this option is not provided, default outputs filename will be `detected_ssrs.bed` and `detected_ssrs_full.bed`. 
Please note that even in the case of no identified SSRs, output files is still created (therefore overwriting any previous file of the same name) but with no content in the file.
The Tool creates 2 files in bed format at the output.
One of its with standard data (more details below), the second (with the suffix `_full.bed`) contains additional information about the found repeats.

**Output format**
The output is a tab-delimited file, with one SSR record per line. The output columns follow the BED format.

|     **Column**     |                                 **Description**                                |
|:------------------:|:------------------------------------------------------------------------------:|
| #Contig            | Chromosome or Sequence Name as specified by the first word in the FASTA header |
| Start              | Start position of SSR in the Chromosome (0-based, included)                    |
| End                | End position of SSR in the Chromosome (0-based, not included)                  |
| Motif              | The motif for which the search was carried out                                 |
| Length             | Total length of identified repeat in bp                                        |
| Strand             | Strand of SSR based on their cyclical variation                                |
| Number of motifs   | Number of times the base motif is repeated                                     |
| Number of inserts* | Number of inserts between motifs in repeat                                     |
| Length of inserts* | Total length of inserts between motifs in repeat                               |
| Sequence           | The found repeat sequence                                                      |

\* -- columns present only in full bed file (*_full.bed).

**`--region or -r`**

**Expects:** _STRING_

**Default:** _None_

This parameter specifies the repeat search region. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire FASTA file will be searched! 0-based format, start is included, end in not included.

**`--cores or -c`**

**Expects:** _INT_

**Default:** _1_

The parameter specifies the number of threads for simultaneous search by contigs and strands. Default is 1 tread.

**`--distance or -d`**

**Expects:** _INT_

**Default:** _0_

The maximum allowable distance between motifs in bp, within which it is assumed that it all falls into the repeat region. Default is 0.

**`--number or -n`**

**Expects:** _INT_

**Default:** _2_

The minimum number of motives in one SSR. Default is 2.

**`--strand or -s`**

**Expects:** _STRING_

**Default:** _both_

This parameter specifies the chain to be searched for. Default is both strand. Use `-s +` to search in a straight chain.
Use `-s -` to search in a reverse chain. Use `-s both` to search in a both chains. 

**`--motif_overlaps`**

**Expects:** _None_

**Default:** _None_

If this option is set it turns on overlapping motif matches.
