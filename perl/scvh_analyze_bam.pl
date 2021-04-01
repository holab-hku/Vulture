#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);

my $samtools = "/home/asdfken/tools/samtools-1.10/samtools";
my $extract_in_cell_gene_barcodes_script = "perl/scvh_extract_in_cell_gene_barcodes.r";

my $output_dir = ".";
my $sample_name = "";
my $threads = 1;

my $STAR_MAPQ_score_for_unique_mapping = 255;

GetOptions(
	't|threads=i' => \$threads,
) or die_usage();

die_usage() if ( (!@ARGV) || (@ARGV > 2));

if (@ARGV == 1) {
	$output_dir = $ARGV[0];
} elsif (@ARGV == 2) {
	$output_dir = $ARGV[0];
	$sample_name = $ARGV[1];
}

sub die_usage {
die "
Usage: scvh_analyze_bam.pl [Options] <output_dir> <sample_name (optional)>

Options:                                                                                                   Defaults
-t/--threads <int>        number of threads to run samtools with                                           [$threads]
";
}

print "Analyze bam in $output_dir\n";

my $date = `date`;
chomp($date);
print "Getting the barcodes of the non-host genes that have potential interaction with host " . "[" . $date . "]\n";
#get which genes need to be included in the in_cell (i.e. EmptyDrops filtered) and result
my %all_in_cell_barcodes;
my %all_in_cell_accessions;
my %in_cell_genes;

my %acc_to_ref_name;

my $sample_name_tag = "";
if ($sample_name ne "") {
	$sample_name_tag = $sample_name . "_";
}

sub get_in_cell_genes {

	my $viral_genes_info = "$output_dir/" . $sample_name_tag . "filtered_matrix_viral_genes_info.txt";
	open my $IN_G, "<$viral_genes_info" or die "can't open $viral_genes_info\n";
	my $header = <$IN_G>;

	my $filtered_matrix = "$output_dir/" . $sample_name_tag . "filtered_matrix.rds";

	while (<$IN_G>) {
		my @line = split("\t", $_);
		my ($gene, $accession, $reference_name) = ($line[0], $line[3], $line[4]);
		chomp($reference_name);
		$acc_to_ref_name{$accession} = $reference_name;

		$in_cell_genes{$gene}{"existing_line"} = chomp($_);
	
		$all_in_cell_accessions{$accession}++;
	}
	close $IN_G;

	my $all_genes = join(",", keys %in_cell_genes);
	my $all_viral_genes_barcodes = `Rscript $extract_in_cell_gene_barcodes_script $filtered_matrix $all_genes`;
	#print "All: $all_viral_genes_barcodes\n";
	my @all_viral_genes_barcodes_ar = split("\n", $all_viral_genes_barcodes);
	for my $viral_genes_barcodes (@all_viral_genes_barcodes_ar) {
		my @viral_gene_barcodes_ar = split(/\s+/, $viral_genes_barcodes);
		my $gene = shift @viral_gene_barcodes_ar;
		print "gene: $gene\n";
		print "bcs: " . join(" ", @viral_gene_barcodes_ar) . "\n\n";
	
		for my $bc (@viral_gene_barcodes_ar) {
			$all_in_cell_barcodes{$bc}++;
			$in_cell_genes{$gene}{"bcs"}{$bc}++;
		}
	}

}

#####################
#####################
get_in_cell_genes();
#####################
#####################
#die;
my $chr_prefix = "chr";

my $STAR_output_dir = $output_dir . "/STARsolo_outs";

unless (-d $STAR_output_dir) {
	die "$STAR_output_dir cannot be found, please run scvh_map_reads.pl first\n";
}

my $analyze_bam_dir = $output_dir . "/analyze_bam";
mkdir $analyze_bam_dir unless (-d $analyze_bam_dir);

$date = `date`;
chomp($date);
print "Looping over STARsolo BAM file to store host read IDS and output non-host BAM " . "[" . $date . "]\n";
#store host (human) read IDs and their relevant info in memory (a hash)
#output non-host reads into a separate BAM
my $BAM = $STAR_output_dir . "/Aligned.sortedByCoord.out.bam";
unless ( (-e $BAM) && (-s $BAM) ) {
	die "$BAM is not present or is empty\n";
}

my $non_host_file = "$analyze_bam_dir/Aligned.sortedByCoord.out.non-host";

my %host_chrs;

my %mapped_to_host_reads;
my %mapped_to_non_host_reads;
my %viral_genes_data;

sub loop_STARsolo_bam {
	open my $OUT_NH, ">$non_host_file.sam" or die "can't open $non_host_file.sam\n";
	#read header lines only first to get the host chrs and to out the non-host chr headers to separate BAM
	open my $IN, "$samtools view -@ $threads -H $BAM |" or die "can't open $BAM\n";
	while (<$IN>) {
		if ($_ =~ /\@HD/) {
			print $OUT_NH $_;
		} elsif ($_ =~ /\@SQ/) {
			my ($ref) = $_ =~ /SN:(\S+)/;
			if ($ref =~ /^$chr_prefix/) {
				$host_chrs{$ref}++;
			} else {
				print $OUT_NH $_;
			}
		} else {
			print $OUT_NH $_;
		}
	}
	close $IN;

	print "host chrs: " . join("\t", keys %host_chrs) . "\n";

	my $mapped_to_host_reads_temp = "$analyze_bam_dir/mapped_to_host_reads.temp";
	open my $OUT_T, ">$mapped_to_host_reads_temp" or die "can't oepn $mapped_to_host_reads_temp\n";
	
	open $IN, "$samtools view -@ $threads $BAM |" or die "can't open $BAM\n";
	while (<$IN>) {
		#line[0] = read, line[2] = ref, line[3] = MAPSTART, line[5] = CIGAR
		my @line = split("\t", $_);
	
		if (exists $host_chrs{$line[2]}) {
			#$mapped_to_host_reads{$line[0]}{"$line[2],$line[3],$line[5]"}++;
			
			#col5 (line[4]) is the MAPQ which according to the STAR manual if the MAPQ is 255, then the mapping is unique and so we don't need to out those
			unless ($line[4] == $STAR_MAPQ_score_for_unique_mapping) {
				print $OUT_T "$line[0]\t$line[2],$line[3],$line[5]\n";
			}
		} else {
			$mapped_to_non_host_reads{$line[0]}{"$line[2],$line[3],$line[5]"}++;
			print $OUT_NH $_;

		}
	}
	close $IN;
	close $OUT_NH;
	close $OUT_T;
	
	system("$samtools view -@ $threads -Sb $non_host_file.sam > $non_host_file.bam");
	system("rm $non_host_file.sam");
	
	print "Total number of non-host read IDs: " . scalar(keys %mapped_to_non_host_reads) . "\n"; 
	
	open my $IN_T, "<$mapped_to_host_reads_temp" or die "can't open $mapped_to_host_reads_temp\n";
	while (<$IN_T>) {
		my @line = split("\t", $_);
		if (exists $mapped_to_non_host_reads{$line[0]}) {
			chomp(@line);
			$mapped_to_host_reads{$line[0]}{$line[1]}++;
		}
	}
	close $IN_T;
	

	system("rm $mapped_to_host_reads_temp");

	print "Total number of host read IDs that have multi-mapping with non-host: " . scalar(keys %mapped_to_host_reads) . "\n";
}

#####################
#####################
loop_STARsolo_bam();
#####################
#####################

$date = `date`;
chomp($date);
print "Looping over non-host BAM file to output reads mapping quality control criteria " . "[" . $date . "]\n";

my $in_cell_reads_out = "$analyze_bam_dir/$sample_name_tag" . "in_cell_reads.tsv";


sub loop_non_host_BAM {
	open my $OUT, ">$in_cell_reads_out" or die "can't open $in_cell_reads_out\n";
	my @out_header = qw(read_ID	CIGAR_string	read_length	UMI	barcode	gene	accession	reference_name 	POS	number_of_times_multi-mapped_to_host	loci_of_mappings_to_host_if_any	number_of_times_multimapped_to_non-host	  	loci_of_mappings_to_non-host);
	print $OUT join("\t", @out_header) . "\n";

	#loop through the non-host BAM and output the final quality control criteria
	#output a non-host.in_cell.bam (only include the non-host SAM headers)
	open my $OUT_NH, ">$non_host_file.in_cell.sam" or die "can't open $non_host_file.in_cell.sam\n";

	open my $IN_B, "$samtools view -@ $threads -h $non_host_file.bam |" or die "can't open $non_host_file.bam\n";
	while (<$IN_B>) {
		if ($_ =~ /^@/) {
			my $acc;
			if ($_ =~ /SN:(\S+)/) {
				$acc = $1;
				if (exists $all_in_cell_accessions{$acc}) {
					print $OUT_NH $_;
				} else {
					next;
				}
			} else {
				next;
			}
		} else {
			#line[0] = read, line[2] = ref, line[3] = MAPSTART, line[5] = CIGAR, line[9] = read_sequence
			my @line = split("\t", $_);
			unless (exists $all_in_cell_accessions{$line[2]}) {
				next;
			}
	
			my ($bc) = $_ =~ /CB:Z:(\S+)/;
			if ($bc) {
				if (exists $all_in_cell_barcodes{$bc}) {
				
					my $read_length = length($line[9]);
					my ($UBZ) = $_ =~ /UB:Z:(\S+)/;
					my ($GNZ) = $_ =~ /GN:Z:(\S+)/;
				
					#read isn't mapped to a gene, will not be counted as in_cell
					if (!$GNZ) {
						next;
					}
				
					print $OUT_NH $_;
					print $OUT "$line[0]\t$line[5]\t$read_length\t$UBZ\t$bc\t$GNZ\t$line[2]\t" . $acc_to_ref_name{$line[2]} . "\t$line[3]\t";
				
					my $number_of_times_multimapped_to_host = 0;
					my @loci_of_mappings_to_host; 
					$loci_of_mappings_to_host[0] = "none";
					if (exists $mapped_to_host_reads{$line[0]}) {
						$number_of_times_multimapped_to_host = scalar (keys %{$mapped_to_host_reads{$line[0]}});
						@loci_of_mappings_to_host = sort keys %{$mapped_to_host_reads{$line[0]}};
					}
					print $OUT "$number_of_times_multimapped_to_host\t" . join(";", @loci_of_mappings_to_host) . "\t";
				
					my $number_of_times_multimapped_to_non_host = scalar (keys %{$mapped_to_non_host_reads{$line[0]}});
					$number_of_times_multimapped_to_non_host -= 1;
					my @loci_of_mappings_to_non_host = sort keys %{$mapped_to_non_host_reads{$line[0]}};
					print $OUT "$number_of_times_multimapped_to_non_host\t" . join(";", @loci_of_mappings_to_non_host) . "\n";
				
				}
			}
		}	
	}
	close $IN_B;

	close $OUT_NH;
	close $OUT;

	system("$samtools view -@ $threads -Sb $non_host_file.in_cell.sam > $non_host_file.in_cell.bam");
	system("$samtools index $non_host_file.in_cell.bam");
	system("rm $non_host_file.in_cell.sam");
}

#####################
#####################
loop_non_host_BAM();
#####################
#####################

#Note: loop_non_host_bam() and loop_STARsolo_bam() require all functions before it to be run, however, output_updated_viral_genes_info() is based totally upon the contents of in_cell_reads_out, and therefore only get_in_cell_genes() needs to be run before output_updated_viral_genes_info()
sub output_updated_viral_genes_info {
	
	my %updated_info;
	
	open my $IN, "<$in_cell_reads_out" or die "can't open $in_cell_reads_out\n";
	my $header_in_cell = <$IN>;
	while (<$IN>) {
		chomp;
		my @line = split("\t", $_);
		my ($read_ID, $CIGAR_string, $read_length, $UMI, $barcode, $gene, $accession, $reference_name, $POS, $number_of_times_multi_mapped_to_host, $loci_of_mappings_to_host_if_any, $number_of_times_multimapped_to_non_host, $loci_of_mappings_to_non_host) = @line;
		
		$updated_info{$accession}{$gene}{"reads"}{$read_ID}++;
		
		#print "$CIGAR_string\n";
		my $num_M=0;
		while ($CIGAR_string =~ /(\d+)M/g) {
			$num_M+=$1;
		}
		#print "$num_M\n";
		my $perc_M = $num_M / $read_length * 100;
		push @{$updated_info{$accession}{$gene}{"perc_M"}}, $perc_M;
		
		push @{$updated_info{$accession}{$gene}{"multi_host"}}, $number_of_times_multi_mapped_to_host;
		
		my @all_loci_of_mappings_to_host_if_any = split(";", $loci_of_mappings_to_host_if_any);
		for my $locus (@all_loci_of_mappings_to_host_if_any) {
			$updated_info{$accession}{$gene}{"loci_multi_host"}{$locus}++;
		}
		
		
		push @{$updated_info{$accession}{$gene}{"multi_non_host"}}, $number_of_times_multimapped_to_non_host;
		
		my @all_loci_of_mappings_to_non_host = split(";", $loci_of_mappings_to_non_host);
		for my $locus (@all_loci_of_mappings_to_non_host) {
			$updated_info{$accession}{$gene}{"loci_multi_non_host"}{$locus}++;
		}
		
		$updated_info{$accession}{$gene}{"POS"}{$POS}++;
		
	}	
	close $IN;
	
	my $updated_viral_genes_info = "$output_dir/" . $sample_name_tag . "filtered_matrix_viral_genes_info.analyze_bam.txt";
	open my $OUT, ">$updated_viral_genes_info" or die "can't open $updated_viral_genes_info\n";
	
	my $viral_genes_info = "$output_dir/" . $sample_name_tag . "filtered_matrix_viral_genes_info.txt";
	open my $IN_G, "<$viral_genes_info" or die "can't open $viral_genes_info\n";
	my $header = <$IN_G>;
	chomp ($header);
	
	my @new_header = qw(number_of_reads	mean_%_M	mean_number_of_times_multi-mapped_to_host	loci_of_mappings_to_host_if_any 	mean_number_of_times_multimapped_to_non-host	loci_of_mappings_to_non-host 	dispersion(number_unique_POS/number_of_reads));
	
	print $OUT $header;
	
	print $OUT "\t" . join("\t", @new_header) . "\n";
	
	while (<$IN_G>) {
		chomp;
		my @line = split("\t", $_);
		my ($gene, $accession, $reference_name) = ($line[0], $line[3], $line[4]);

		print $OUT $_;
		
		if (exists $updated_info{$accession}) {
			if (exists $updated_info{$accession}{$gene}) {
				#number_of_reads
				print $OUT "\t" . scalar(keys %{$updated_info{$accession}{$gene}{"reads"}});
				
				#mean_%_M
				my $mean_perc_M = mean(@{$updated_info{$accession}{$gene}{"perc_M"}});
				print $OUT "\t" . $mean_perc_M;
				
				#mean_number_of_times_multi-mapped_to_host
				my $mean_multi_host = mean(@{$updated_info{$accession}{$gene}{"multi_host"}});
				print $OUT "\t" . $mean_multi_host;
				
				#loci_of_mappings_to_host_if_any [locus (num reads);]
				my @loci_host;
				for my $locus (keys %{$updated_info{$accession}{$gene}{"loci_multi_host"}}) {
					my $locus_num_reads = "$locus (" . $updated_info{$accession}{$gene}{"loci_multi_host"}{$locus} . ")";
					push @loci_host, $locus_num_reads;
				}
				print $OUT "\t" . join(";", @loci_host);
				
				#mean_number_of_times_multi-mapped_to_non-host
				my $mean_multi_non_host = mean(@{$updated_info{$accession}{$gene}{"multi_non_host"}});
				print $OUT "\t" . $mean_multi_non_host;
				
				#loci_of_mappings_to_non-host [locus (num reads);]
				my @loci_non_host;
				for my $locus (keys %{$updated_info{$accession}{$gene}{"loci_multi_non_host"}}) {
					my $locus_num_reads = "$locus (" . $updated_info{$accession}{$gene}{"loci_multi_non_host"}{$locus} . ")";
					push @loci_non_host, $locus_num_reads;
				}
				print $OUT "\t" . join(";", @loci_non_host);
				
				my $dispersion = scalar(keys %{$updated_info{$accession}{$gene}{"POS"}}) / scalar(keys %{$updated_info{$accession}{$gene}{"reads"}});
				print $OUT "\t" . $dispersion . "\n"; 
				
			} else {
				die "$gene not in updated info\n";
			}
		} else {
			die "$accession not in updated info\n";
		}
		
		#die;
	}
	close $IN_G;
	close $OUT;
}

output_updated_viral_genes_info();


#(gene	UMI_counts	median_num_human_genes_expressed	accession	reference_name)<- existing	total no. of reads	mean_%_M	mean no. times multi-mapped to host		mean no. times multi-mapped to non-host	loci of maps to host (if any)	loci of maps to non-host		number_unique_MAPSTART / number_of_reads
sub mean {
	my @ar = @_;
	my $tot;
	for my $ele (@ar) {
		$tot += $ele;
	}
	my $mean = $tot / scalar(@ar);
	return $mean;
}

