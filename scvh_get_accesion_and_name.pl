#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);

my $output_dir = ".";
my $genome_dir = "\"vh_genome_dir\"";

GetOptions(
	'o|output-dir=s' => \$output_dir

) or die_usage();

die_usage() unless @ARGV == 1;

#my ($genome_dir, $R2, $R1) = @ARGV;
$genome_dir = $ARGV[0];

sub die_usage {
die "
Usage: scvh_map_reads.pl [Options] <vh_genome_dir> <R2> <R1>

Options:                                                                                                   Defaults
-o/--output-dir	<string>  the output directory                                                             [./]
";
}

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
# index_STAR_genome_if_nec();
# run_STAR();

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