############################################################
#
# RS_aligner.pl
# - Reap-Seq Aligner: An aligner to align antibody oligo reads to a dictionary
#    
# by Kelvin Zhang
# 
# (c) 2017 Kelvin Zhang
# Version 06.06.2017
# 
############################################################

# HOW TO RUN
# perl RS_aligner.pl YOUR_SAM_FILE OUTPUT_PREFIX 
# Example: 
# perl RS_aligner.pl -i test.sam -d Antibody_oligo_sequence_dictionary.fa -oname test

#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw<first>;
use List::Util qw(reduce);
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(uniq);

# version
my $VERSION = "1.0";
print "RS_aligner version $VERSION\n";

# get options
my ($INPUT_FILE, $OPREFIX, $DICT, $DISTANCE);
 
GetOptions(
	'INPUT_FILE|i=s'  => \$INPUT_FILE,
	'OPREFIX|oname=s'  => \$OPREFIX,
	'DICT|d=s'  => \$DICT,
	'DISTANCE|m=s'  => \$DISTANCE,
);


my $numArgs = $#ARGV + 1;

my @arguments;
foreach my $argnum (0 .. $#ARGV){
     push(@arguments, $ARGV[$argnum]);
}
 
# input
my @dictionary = read_file($DICT);
open(INFILE2, "<$INPUT_FILE");
 
#output
my $outputfile = "$OPREFIX"."_matrix.txt"; # a file used to generate a matrix in which each row is oligo and each column is cell
open(OUT, ">$outputfile");

my $outputfile2 = "$OPREFIX"."_alignment_all.txt"; # details about mapping of each oligo
open(OUT2, ">$outputfile2");

my $outputfile3 = "$OPREFIX"."_counts_by_cellBarcode.txt"; # read counts of each cell
open(OUT3, ">$outputfile3");

# reference hash
my %reference_hash;
for(my $i=0; $i < @dictionary; $i++){
	my $key = $dictionary[$i];
#	print "key is: $key\n";
	$key =~ s/>//;
	$reference_hash{$key} = $dictionary[$i+1];
	$i = $i + 1;
}

 
my %countTable;
my %cb_hash_total_read_count;

print OUT2 "oligo_mapped\tread\tcell_barcode\tUMI\n";
print OUT3 "cell_barcode\ttotal_count\n";

while (<INFILE2>) {
	chomp;
	if($_ =~ /CB:Z/ && $_ =~ /UB:Z/){
		my @items = split("\t", $_);
#		print "$items[9]\n"; # read sequence
#		my $read = substr $items[9], 0, 8; # first 8 bp of read
		my $read = $items[9]; # first 8 bp of read
		my $index = align_reference($read, $DISTANCE);
#		print "index is: $index\n";
		my $cb;
		my $umi;
		for(my $j=0; $j < @items; $j++){
			if($items[$j] =~ /^CB/){
				$cb = $items[$j];
				$cb =~ s/CB:Z://;
#				print "$cb\n";
			}
			if($items[$j] =~ /^UB/){
				$umi = $items[$j];
				$umi =~ s/UB:Z://;
#				print "$umi\n";
			}
		}
		
		$cb_hash_total_read_count{$cb}++;
		 
		print OUT2 "$index\t$read\t$cb\t$umi\n";
		 
		$countTable{$index}{$cb}{$umi} = 1;
		
	}
}
	
print "The processing of sam file is done.\n";	

foreach my $cb (reverse sort { $cb_hash_total_read_count{$a} <=> $cb_hash_total_read_count{$b} } keys %cb_hash_total_read_count) {
	print OUT3 "$cb\t$cb_hash_total_read_count{$cb}\n";
#	print "$cb\t$cb_hash_total_read_count{$cb}\n";
}

#print Dumper \%matrix;

my %matrix;
my @cell_list;
my @gene_list;

foreach my $gene (keys %countTable){ # by gene
	foreach my $cell (keys %{$countTable{$gene} }) {
		my $count = keys %{$countTable{$gene}{$cell}};
		my @list = keys %{$countTable{$gene}{$cell}}; 
		my $value = "$gene\t$cell\t$count";
#		print OUT "$gene\t$cell\t$count\n";
		 
		if($gene ne "Not_mapped"){
			push(@cell_list, $cell);
			push(@gene_list, $gene);
			$matrix{$gene}{$cell} = $count;
		}
	}
}
	
my @cell_list_new = uniq @cell_list;
my @gene_list_new = uniq @gene_list;

my $count_cell = @cell_list_new;
my $count_gene = @gene_list_new;

print "The number of detected cells are: $count_cell\n";
print "The number of detected proteins are: $count_gene\n";

my $header = join("\t", @cell_list_new);
print OUT "Gene\t$header\n";

for(my $i=0; $i < @gene_list_new; $i++){
	my @list;
	for(my $j=0; $j < @cell_list_new; $j++){
		if($matrix{$gene_list_new[$i]}{$cell_list_new[$j]}){
			push(@list, $matrix{$gene_list_new[$i]}{$cell_list_new[$j]});
		}else{
			push(@list, 0);
		}
	}
	my $string = join("\t", @list);
#	print "$gene_list_new[$i]\t$string\n";
	print OUT "$gene_list_new[$i]\t$string\n";
}


exit;

# align read to reference hash
sub align_reference {
	my ($seq, $mismatch) = @_;
	
	my $index = "Not_mapped";
	
	foreach my $key (keys %reference_hash){
		my $reference_seq = $reference_hash{$key};
		
#		print "query is: $seq\n";
#		print "reference is: $reference_seq\n";
		
		if(mismatch($seq, $reference_seq, $mismatch)){ # use hamming distance to allow one mismatch
#		if($seq eq $reference_seq){ # no mismatch
			$index = $key;
		}
	}
	
	return $index;
}
	
# function to find match max_distance = mismatch
sub mismatch {
	my ($seq1, $seq2, $max_distance) = @_;
	
	if($max_distance > 1){
		for my $offset (0..length($seq2)-length($seq1)) {
			my $hd = hd($seq1,substr($seq2,$offset,length($seq1)));
			if ($hd <= $max_distance) {
				return 1;
			}else{
				return 0;
			}
		}
	}else{
		if($seq1 eq $seq2){
			return 1;
		}else{
			return 0;
		}
	}
}
	
# assumes byte mode
sub hd {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

# read file
sub read_file {
	my ($filename) = @_;
	
	open(INPUTFILE1, $filename) or die "can not open the file!\n";
	my @input = <INPUTFILE1>;
	chomp(@input);
	return @input;
}
