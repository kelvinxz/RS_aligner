# RS_aligner
An aligner to align antibody oligo reads to an antibody dictionary

## Using RS_alinger

Usage: perl RS_aligner.pl -i YOUR_SAM_FILE -d DICTIONARY -oname OUTPUT_PREFIX -m MISMATCH

## Arguments:
* **-i**      A sam file processed by the 10X CellRanger pipeline.
* **-d**      An antibody dictionary (fasta format)
* **-oname**  A output prefix.
* **-m**     Numbers of mismatches allowed during alignment.

## An example: 

perl RS_aligner.pl -i MyTest.sam -d Antibody_oligo_sequence_dictionary.fa -oname MyTest -m 1
