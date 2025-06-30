#  Deduplicate TE library 

## Introduction
The script is written to remove overlapping TE regions in the EDTA results to fit the buildsummary.pl statistical script. More details see https://github.com/Dfam-consortium/RepeatMasker/issues/343

The criteria for deletion prioritize the retention of intact LTRs (if --intact-gff present), followed by TEs identified by homology, then Unknown TEs identified by homology, next TEs identified by structural methods, and finally Unspecified Repeats and Tandem Repeats identified by TRF.

If overlapping TEs have the same priority level, the longer one will be retained preferentially.

## Usage
```bash
# Duplicate TE
python TEdedup.py INBED OUTBED --intact-gff Intact.gff3

# Generate the .out file required for the buildsummary.pl script.
perl -nle '
    use feature "declared_refs";
    no warnings "experimental::declared_refs";
    use feature "refaliasing";
    no warnings "experimental::refaliasing";

    my ($chr, $s, $e, $anno, $dir, $supfam) = (split)[0,1,2,3,8,12];
    print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"
    ' OUTBED > OUT

# Summary
perl buildSummary.pl -maxDiv 40 -genome_size GenomeSize -seq_count SeqNum OUT > SUM 2>/dev/null
```

## Example
```bash
python TEdedup.py test_data/reclassify.chr1.bed test_data/reclassify.chr1.split.bed --intact-gff test_data/intact.chr1.gff3

perl -nle '
    use feature "declared_refs";
    no warnings "experimental::declared_refs";
    use feature "refaliasing";
    no warnings "experimental::refaliasing";

    my ($chr, $s, $e, $anno, $dir, $supfam) = (split)[0,1,2,3,8,12];
    print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"
    ' test_data/reclassify.chr1.split.bed > test_data/reclassify.chr1.split.out

perl utils/buildSummary.pl -maxDiv 40 -genome_size 58642971 -seq_count 1 test_data/reclassify.chr1.split.out > test_data/reclassify.chr1.split.summary 2>/dev/null
```