#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);

my $output_dir = ".";
my $threads = 1;
my $soloStrand = "Reverse"; #Reverse is used for 10x 5' protocol, while Forward is used for 10x 3' protocol
my $genome_dir = "\"vh_genome_dir\"";
my $barcodes_whitelist = "$genome_dir/737K-august-2016.txt";

GetOptions(
	'o|output-dir=s' => \$output_dir,
	't|threads=i' => \$threads,
	's|soloStrand=s' => \$soloStrand,
	'w|whitelist=s' => \$barcodes_whitelist,
) or die_usage();

die_usage() unless @ARGV == 3;

#my ($genome_dir, $R2, $R1) = @ARGV;
$genome_dir = $ARGV[0];
my ($R2, $R1) = ($ARGV[1], $ARGV[2]);

sub die_usage {
die "
Usage: scvh_map_reads.pl [Options] <vh_genome_dir> <R2> <R1>

Options:                                                                                                   Defaults
-o/--output-dir	<string>  the output directory                                                             [./]
-t/--threads <int>        number of threads to run STARsolo with                                           [$threads]
-s/--soloStrand <string>  STARsolo param: Reverse or Forward used for 10x 5' or 3' protocol, respectively  [$soloStrand]
-w/--whitelist <string>   STARsolo param --soloCBwhitelist                                                 [<$barcodes_whitelist>]
";
}

if ($barcodes_whitelist ne "\"vh_genome_dir\"/737K-august-2016.txt") {
	
} else {
	$barcodes_whitelist = "$genome_dir/737K-august-2016.txt";
}

#params to be removed
#definitely needed, can use static version for now, but eventually user should pre-install STAR
my $STAR = "/content/gdrive/MyDrive/STAR-2.7.8a/bin/Linux_x86_64_static/STAR";

#not necessary unless analyze_BAM() is called, which currently it isn't
my $samtools = "/home/asdfken/tools/samtools-1.10/samtools";

#params not planned to be user options (for now)
my $host_species = "human";
my $host_ref_genome = "hg38";
my $chr_prefix = "chr";
my $virus_database = "viruSITE";
my $use_removed_amb_viral_exon = "T";
my $removed_amb_viral_exon_tag = "removed_amb_viral_exon";
my @STAR_index_files = qw(SAindex SA Genome genomeParameters.txt chrStart.txt chrNameLength.txt chrName.txt chrLength.txt); #as of v2.7.5a
my $soloUMIfiltering = "-";
my $max_multi = 20;
my $soloCellFilter = "None";
my $soloBarcodeReadLength = 0;

mkdir $output_dir unless (-d $output_dir);

my %human_chrs;
my %acc_to_name;

my ($fa, $gtf) = check_genome_fasta_gtf_present();
get_reference_names_and_accessions();
out_gene_to_accession_and_name();
index_STAR_genome_if_nec();
run_STAR();

sub out_gene_to_accession_and_name {
	my $intermediate_output_dir = "$output_dir/intermediate_files";
	mkdir $intermediate_output_dir unless (-d $intermediate_output_dir);
	my $gene_to_accession_and_name_out_file = $intermediate_output_dir . "/gene_to_accession_and_name.txt";
	open my $OUT, ">$gene_to_accession_and_name_out_file" or die "can't open $gene_to_accession_and_name_out_file\n";
	
	print $OUT "gene\t" . "accession\t" . "reference_name\n";
	
	open my $IN, "<$gtf" or die "can't open $gtf\n";
	
	my %already_out;
	
	while (<$IN>) {
		my @line = split("\t", $_);
		my ($ref, $feature) = ($line[0], $line[2]);
		if (exists $human_chrs{$ref}) {
			next;
		} elsif (exists $acc_to_name{$ref}) {
			unless ($feature eq "exon") {
				next;
			}
			
			my $name = $acc_to_name{$ref};
			my ($gene_name) = $_ =~ /gene_name \"(.*?)\";/;
			my ($gene_id) = $_ =~ /gene_id \"(.*?)\";/;
			my $final_name;
			if ($gene_id && $gene_name) {
				$final_name = $gene_name;
			} elsif ($gene_id) {
				$final_name = $gene_id;
			} else {
				print "Error, can't find gene_name or gene_id for $_";
			}
			
			if (exists $already_out{$final_name}) {
				next;
			}
			
			print $OUT $final_name . "\t" . $ref . "\t" . $name . "\n";
			
			$already_out{$final_name}++;
			
		}
	}
	close $IN;
	close $OUT;
	
}

sub analyze_BAM {
	my $STAR_output_dir = shift;
	my $BAM = $STAR_output_dir . "/Aligned.sortedByCoord.out.bam";
	unless ( (-e $BAM) && (-s $BAM) ) {
		die "$BAM is not present or is empty\n";
	}
	
	my %mapped_to_human_reads;
	my %viral_genes_data;
	
	open my $IN, "$samtools view -@ $threads $BAM |" or die "can't open $BAM\n";
	while (<$IN>) {
		#line[0] = read, line[2] = ref, 
		my @line = split("\t", $_);
		
		if (exists $human_chrs{$line[2]}) {
			$mapped_to_human_reads{$line[0]}++;
		} else {
			my ($GXZ) = $_ =~ /GX:Z:(\S+)/;
			my ($GNZ) = $_ =~ /GN:Z:(\S+)/;
			my ($UBZ) = $_ =~ /UB:Z:(\S+)/;
			
			if ($UBZ && $GXZ) {
				$viral_genes_data{$GXZ}{"GNZ"} = $GNZ;
				$viral_genes_data{$GXZ}{"UBZ"}{$UBZ}++;
				$viral_genes_data{$GXZ}{"acc"} = $line[2];
				
				
				if (exists $mapped_to_human_reads{$line[0]}) {
					$viral_genes_data{$GXZ}{"mth_reads"}++;
				}
				
				$viral_genes_data{$GXZ}{"tot_reads"}++;
			}
		}
	}
	close $IN;
	
	my $output_file = "$output_dir/mapped_viral_genes.txt";
	open my $OUT, ">$output_file" or die "can't open $output_file\n";
	my @out_cols = qw(gene_id gene_name accession reference_name UMI_count_unfiltered read_count percent_reads_multimapped_to_host);
	print $OUT join("\t", @out_cols) . "\n"; 
	for my $gene_id (sort keys %viral_genes_data) {
		my $acc = $viral_genes_data{$gene_id}{"acc"};
		my $name = $acc_to_name{$acc};
		print $OUT $gene_id . "\t" . $viral_genes_data{$gene_id}{"GNZ"} . "\t"  . $acc . "\t" . $name . "\t";
		print $OUT scalar(keys %{$viral_genes_data{$gene_id}{"UBZ"}}) . "\t" . $viral_genes_data{$gene_id}{"tot_reads"} . "\t";
		my $percent_reads_multimapped_to_host;
		if (!exists $viral_genes_data{$gene_id}{"mth_reads"}) {
			$percent_reads_multimapped_to_host = 0;
		} else {
			$percent_reads_multimapped_to_host = $viral_genes_data{$gene_id}{"mth_reads"} / $viral_genes_data{$gene_id}{"tot_reads"} * 100;
		}
		print $OUT $percent_reads_multimapped_to_host . "\n";
		
	}
	close $OUT;
}

sub run_STAR {
	
	for my $R ($R2, $R1) {
		unless ( (-e $R) && (-s $R) ) {
			die "$R is not present or is empty\n";
		}
	}
	
	my $readFilesCommand;
	
	if ( ($R2 =~ /.*\.gz$/) && ($R1 =~ /.*\.gz$/) ) {
		$readFilesCommand = "zcat";
	} else {
		print ".gz not detected in file suffix of both R1 and R2, will assume both are uncompressed.\n";
		$readFilesCommand = "-";
	}
	
	my $STAR_output_dir = $output_dir . "/STARsolo_outs/";
	mkdir $STAR_output_dir unless (-d $STAR_output_dir);
	
	my $run_STAR = "$STAR --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";

	system("$run_STAR");
	
	#analyze_BAM($STAR_output_dir);
}

sub index_STAR_genome_if_nec {
	
	my $index_STAR = "F";
	
	for my $STAR_index_file (@STAR_index_files) {
		if (-e "$genome_dir/$STAR_index_file") {
			if (-s "$genome_dir/$STAR_index_file") {
			} else {
				print "$genome_dir/$STAR_index_file exists but is empty, will index STAR genome\n";
				$index_STAR = "T";
				last;
			}
		} else {
			print "Can't find $genome_dir/$STAR_index_file, will index STAR genome\n";
			$index_STAR = "T";
			last;
		}
	}
	
	if ($index_STAR eq "F") {
		print "STAR genome index files all found, will proceed with the existing index\n";
	} elsif ($index_STAR eq "T") {
		my $generate_genome = "$STAR --runThreadN $threads --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $fa";
		system("$generate_genome");
	}
}

sub get_reference_names_and_accessions {

	print "Getting reference names and accessions...\n";
	open my $IN, "cat $fa | grep \">\" |" or die "can't open $fa\n";
	while (<$IN>) {
		if ($_ =~ /^>($chr_prefix\S+)/) {
			$human_chrs{$1}++;
		} else {
			if ($_ =~ /^>(\S+)\s+(.*)/) {
				my $acc = $1;
				my $name = $2;
				chomp($name);
				$acc_to_name{$acc} = $name;
			}
		}
	}
	close $IN;
}

sub check_genome_fasta_gtf_present {
	my $genome_fa = $genome_dir . "/$host_species" . "_host_viruses.$virus_database.with_" . $host_ref_genome . ".fa";
	unless ( (-e $genome_fa) && (-s $genome_fa) ) {
		die "$genome_fa is not present or is empty\n";
	}
	
	my $genome_gtf = $genome_dir . "/$host_species" . "_host_viruses.$virus_database.with_" . $host_ref_genome;
	if ($removed_amb_viral_exon_tag) {
		$genome_gtf .= ".$removed_amb_viral_exon_tag.gtf";
	} else {
		$genome_gtf .= ".gtf";
	}
	unless ( (-e $genome_gtf) && (-s $genome_gtf) ) {
		die "$genome_gtf is not present or is empty\n";
	}
	
	return ($genome_fa, $genome_gtf);
}