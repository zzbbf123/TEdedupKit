#  Deduplicate TE library 

## Introduction
Deduplicate overlap TE in TE library

## Usage
```bash
python TEdedup.py INBED OUTBED --intac-gff Intact.gff3
```

## Example
```bash
python TEdedup.py test_data/reclassify.chr1.bed test_data/reclassify.chr1.split.bed --intac-gff test_data/intact.chr1.gff3
```