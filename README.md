# Bivartect

**Accurate and memory-saving detection of genomic structural variations via direct comparison of sequence reads**

Last updated: 2018-08-02

We present Bivartect, a genomic structural variant caller that directly compares sequence reads generated by high-throughput sequencing. Bivartect achieves memory saving by keeping only a small part of the suffixes of input reads in memory. Using simulated benchmark data and real genome editing data, Bivartect was at least comparable to state-of-the-art algorithms in predictive performance, especially better for single nucleotide variant detection.

## Installation
* Bivartect (ver. 1.1.9) (**bivartect-1.1.9.tar.gz**) in C++ program

### Requirements
* C++11 or later

### Install on Linux and macOS
Type the followings in your terminal:
```
$ tar zxf bivartect-1.1.9.tar.gz
$ cd bivartect-1.1.9
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
 -d     Filtering depth (int 10...32 [22])
 -c     Read count cutoff.
        In a breakpoint cluster, 
        IF max(predictedNormalReadCount, predictedMutatedReadCount) < c 
        THEN omit the breakpoint because of low quality. (int 1...100 [8])
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

## Data
* Simulated benchmark data

## Reference
Keisuke Shimmura, Yuki Kato and Yukio Kawahara,
**Bivartect: accurate and memory-saving detection of genomic structural variations via direct comparison of sequence reads**.
*in preparation*.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/)  
*Graduate School of Medicine, Osaka University, Japan*
