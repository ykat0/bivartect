# Bivartect

[![DOI](https://zenodo.org/badge/143263549.svg)](https://zenodo.org/badge/latestdoi/143263549)

### Accurate and memory-saving breakpoint detection by direct read comparison

Last updated: 2020-01-23

We present Bivartect, a genomic structural variant caller that directly compares sequence reads generated by high-throughput sequencing. Bivartect achieves memory saving by keeping only a small part of the suffixes of input reads in memory. Using simulated benchmark data and real genome editing data, Bivartect outperformed the state-of-the-art small variant callers in low false positive detection of single nucleotide variants.

## Installation
* Bivartect (ver. 1.1.10) (**bivartect-1.1.10.tar.gz**) in C++ program

### Requirements
* C++11 or later

### Install on Linux and macOS
Type the followings in your terminal:
```
$ tar zxf bivartect-1.1.10.tar.gz
$ cd bivartect-1.1.10
$ ./configure
```
or  
```
$ ./configure CXXFLAGS='-std=c++11 -pthread'
```
If you would like to install your local directory,

```
$ ./configure --prefix=/path/to/local_dir
```
Then,

```
$ make
$ sudo make install
```

## Usage
```
For single-end reads:
$ bivartect -3 <normal.fastq> <tumor.fastq> <output.fastq>

For paired-end reads:
$ bivartect -5 <normal_1.fastq> <normal_2.fastq> <tumor_1.fastq> <tumor_2.fastq> <output.fastq>

General options:
 -n     Path to the normal FASTQ (string [necessary])
 -N     Path to the normal reversed FASTQ (string)
 -m     Path to the mutated FASTQ (string [necessary])
 -M     Path to the mutated reversed FASTQ (string)
 -o     Path to the output FASTQ (string)
 -a     Output multi-FASTA instead of FASTQ (bool [false])
 -s     Input FASTQ is strand-specific (bool [false])
 -d     Filtering depth (int 10...32 [24])
 -c     Read count cutoff.
        In a breakpoint cluster, 
        IF max(predictedNormalReadCount, predictedMutatedReadCount) < c 
        THEN omit the breakpoint because of low quality. (int 1...100 [6])
 -x     Analysis division rate (int 1,4,16,64...1024 [64])
 -t     Using thread count. Set 0 to use hardware maximum threads (int 0... [0])
 -r     Path to the output detail overview text file (string)

Alias options:
 -2     = -n -m
 -3     = -n -m -o
 -4     = -n -N -m -M
 -5     = -n -N -m -M -o

Examples:
$ bivartect -x 16 -d 30 -c 6 -n <normal.fastq> -m <tumor.fastq> -o <output.fastq>
$ bivartect -3 <normal.fastq> <tumor.fastq> <output.fastq> -c 4
$ bivartect -5 <normal_1.fastq> <normal_2.fastq> <tumor_1.fastq> <tumor_2.fastq> <output.fastq>
$ bivartect -2 <normal.fastq> <tumor.fastq> -r <output.txt>
```
## Pipeline
The standard use of Bivartect is illustrated with the following steps:

### Step 1: run Bivartect to get consensus normal FASTQ reads whose mutated counterparts are predicted to have breakpoints
```
$ bivartect -5 <normal_1.fastq> <normal_2.fastq> <tumor_1.fastq> <tumor_2.fastq> <out.fastq>
```

### Step 2: map FASTQ reads onto a reference genome with BWA-backtrack
```
$ bwa aln <index_prefix> <out.fastq> > <out.sai>
$ bwa samse -f <out.sam> <index_prefix> <out.sai> <out.fastq>
```

### Step 3: convert SAM alignments into predicted VCF variants with their genomic locations
```
$ ./sam2vcf.py <out.sam> <reference.fa.gz> > <out.vcf> 
```

## Data
* Simulated benchmark FASTQ data used in this work are available [HERE](http://www.med.osaka-u.ac.jp/pub/rna/ykato/project/bivartect/).

## Reference
Keisuke Shimmura, Yuki Kato and Yukio Kawahara,
**Bivartect: accurate and memory-saving breakpoint detection by direct read comparison**,
*Bioinformatics*, in press.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/)  
*Graduate School of Medicine, Osaka University, Japan*
