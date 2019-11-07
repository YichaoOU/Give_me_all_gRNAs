## Introduction

There are many online or offline CRISPR programs that can give  you a list of gRNAs given input regions. However, most likely they will filter out low-confidence gRNAs. For a very small region, such as 1kb to 5kb, you might want to use all gRNAs in this region. In this case, available tools are limited. Cas-OFFinder can do this, yet you need a little bit more steps than just running it. 

This small program makes it easier to use Cas-OFFinder and output all gRNAs in your input regions.

## Input

Any number of columns are OK. The first 3 columns have to be chr, start, end

```
chr6	135376182	135376982
chr6	135375287	135375739
```

## Output

1. candidate_gRNA.bed

This is a bed-6 format. The columns are chr, start, end, gRNA_sequence, GC%, strand.

2. candidate_gRNA.bed.off_targets.info.csv

In addition to bed-6 file above, this file provides two more columns, namely Number of off targets (exact match) and off targets coordinates.

## Installation

You need to have the following program installed.

1. bedtools
2. cas-offinder
3. pandas - Python library

## Usage

`python find_all_gRNA.py -f input.bed`
