#  Deduplicate TE library 

## Introduction
The script is written to remove overlapping TE regions in the EDTA results to fit the buildsummary.pl statistical script. More details see https://github.com/Dfam-consortium/RepeatMasker/issues/343

The criteria for deletion prioritize the retention of intact LTRs (if --intact-gff present), followed by TEs identified by homology, then Unknown TEs identified by homology, next TEs identified by structural methods, and finally Unspecified Repeats and Tandem Repeats identified by TRF.

If overlapping TEs have the same priority level, the longer one will be retained preferentially.

## Usage
```bash
python TEdedup.py INBED OUTBED --intac-gff Intact.gff3
```

## Example
```bash
python TEdedup.py test_data/reclassify.chr1.bed test_data/reclassify.chr1.split.bed --intac-gff test_data/intact.chr1.gff3
```