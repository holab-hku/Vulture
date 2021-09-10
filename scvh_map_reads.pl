#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);
use Cwd 'abs_path';
my $dir = abs_path($0);
my @dir_str = split('/', "$dir");
my $code_dir = join "/", @dir_str[0 .. $#dir_str-1];
my $output_dir = ".";
my $threads = 1;
my $soloStrand = "Reverse"; #Reverse is used for 10x 5' protocol, while Forward is used for 10x 3' protocol
my $genome_dir = "\"vmh_genome_dir\"";
my $barcodes_whitelist = "$genome_dir/737K-august-2016.txt";
#definitely needed, can use static version for now, but eventually user should pre-install STAR
my $ram = 8;
my $alignment = "STAR";
my $technology = "10XV2";
my $virus_database = "viruSITE.NCBIprokaryotes";
my $pseudoBAM = "";
my $soloMultiMappers  = "EM";
my $soloFeatures = "Gene";
my $outSAMtype = "BAM SortedByCoordinate";
my $soloCBstart = 1;
my $soloCBlen = 16;
my $soloUMIstart = 17;
my $soloUMIlen = 10;
my $EXE = "";

GetOptions(
	'o|output-dir=s' => \$output_dir,
	't|threads=i' => \$threads,
	'r|ram=i' => \$ram,
	'd|database=s' => \$virus_database,
	'e|exe=s' => \$EXE,
	's|soloStrand=s' => \$soloStrand,
	'w|whitelist=s' => \$barcodes_whitelist,
	'a|alignment=s' => \$alignment,
	'x|technology=s' => \$technology,
	'f|soloFeature=s' => \$soloFeatures, 
	'ot|outSAMtype=s' => \$outSAMtype,
	'mm|soloMultiMappers=s' => \$soloMultiMappers,
	'pseudoBAM' => \$pseudoBAM,
	'soloCBstart=i' => \$soloCBstart,
	'soloCBlen=i' => \$soloCBlen,
	'soloUMIstart=i' => \$soloUMIstart,
	'soloUMIlen=i' => \$soloUMIlen

) or die_usage();

die_usage() unless @ARGV == 3;

#my ($genome_dir, $R2, $R1) = @ARGV;
$genome_dir = $ARGV[0];
my ($R2, $R1) = ($ARGV[1], $ARGV[2]);


if($EXE eq "") {
	if($alignment eq "STAR"){
		$EXE ="STAR"; 
	}elsif($alignment eq "KB"){
		$EXE ="kb"; 
	}elsif($alignment eq "Alevin"){
		$EXE ="salmon"; 
	}elsif($alignment eq "CellRanger"){
		$EXE ="cellranger"; 
	}
}

sub die_usage {
die "
Usage: scvh_map_reads.pl [Options] <vmh_genome_dir> <R2> <R1>

Options:                                                                                                                                Defaults
-o/--output-dir	<string>   the output directory                                                                                          [./]
-t/--threads <int>         number of threads to run alignment with                                                                       [<$threads>]
-d/--database <string>     select virus or virus and prokaryotes database, can be 'viruSITE' or 'viruSITE.NCBIprokaryotes'               [<$virus_database>]
-e/--exe <string>          executable command or stand alone executable path of the alignment tool                                       [<$EXE>]
-s/--soloStrand <string>   STARsolo param: Reverse or Forward used for 10x 5' or 3' protocol, respectively                               [<$soloStrand>]
-w/--whitelist <string>    STARsolo param --soloCBwhitelist                                                                              [<$barcodes_whitelist>]
-r/--ram <int>             limitation of RAM usage. For STARsolo, param: limitGenomeGenerateRAM unit by GB                               [<$ram>]
-f/--soloFeatures <string> STARsolo param:  See --soloFeatures in STARsolo manual                                                        [<$soloFeatures>]
-ot/--outSAMtype <string>  STARsolo param:  See --outSAMtype in STARsolo manual                                                          [<$outSAMtype>]
-mm/--soloMultiMappers <string>  STARsolo param:  See --soloMultiMappers in STARsolo manual                                              [<$soloMultiMappers>]
-a/--alignment <string>    Select alignment methods: 'STAR', 'KB', 'Alevin', or 'CellRanger'                                             [<$alignment>]
-v/--technology <string>   KB param:  Single-cell technology used (`kb --list` to view)                                                  [<$technology>]
--soloCBstart <string>  STARsolo param:  See --soloCBstart in STARsolo manual                                                            [<$soloCBstart>]
--soloCBlen <string>  STARsolo param:  See --soloCBlen in STARsolo manual                                                                [<$soloCBlen>]
--soloUMIstart <string>  STARsolo param:  See --soloUMIstart in STARsolo manual                                                          [<$soloUMIstart>]
--soloUMIlen <string>  STARsolo param:  See --soloUMIlen in STARsolo manual                                                              [<$soloUMIlen>]
";
}

if ($barcodes_whitelist ne "\"vmh_genome_dir\"/737K-august-2016.txt") {
	
} else {
	$barcodes_whitelist = "$genome_dir/737K-august-2016.txt";
}

#params to be removed
#not necessary unless analyze_BAM() is called, which currently it isn't
my $samtools = "/home/asdfken/tools/samtools-1.10/samtools";

if ($alignment eq "STAR") {
	$ram = $ram * 1073741274;
}

print $virus_database;


#params not planned to be user options (for now)
my $host_species = "human";
my $host_ref_genome = "hg38";
my $chr_prefix = "chr";
my $use_removed_amb_viral_exon = "T";
my $removed_amb_viral_exon_tag = "removed_amb_viral_exon";
my @STAR_index_files = qw(SAindex SA Genome genomeParameters.txt chrStart.txt chrNameLength.txt chrName.txt chrLength.txt); #as of v2.7.5a
my @KB_index_files = qw(transcriptome.idx cdna.fa transcripts_to_genes.txt); #as of v2.7.5a

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

if ($alignment eq "STAR") {
	index_STAR_genome_if_nec();
	run_STAR();
}elsif ($alignment eq "KB") {
	KBref_if_nec();
	run_KB();
	convert_h5ad_to_10x();
	#convert_BUS_to_text();
}elsif ($alignment eq "Alevin") {
	Alevin_if_nec();
	run_Alevin();

}elsif ($alignment eq "CellRanger") {
	CRref_if_nec();
	run_CR();

}


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
	
	my $STAR_output_dir = $output_dir . "/alignment_outs/";
	mkdir $STAR_output_dir unless (-d $STAR_output_dir);

	my $outSAMattributes = "CR CY CB UR UY UB sM GX GN NH";

	my @strSAMtype = split(" ", $outSAMtype);
	
	if( 'Unsorted' ~~ @strSAMtype ){
		$outSAMattributes = "CR CY UR UY sM GX GN NH";
	}
	#my $run_STAR = "$EXE --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";
	my $run_STAR = "$EXE --runThreadN $threads --genomeDir $genome_dir --limitGenomeGenerateRAM $ram --outSAMtype $outSAMtype --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloFeatures $soloFeatures --soloMultiMappers $soloMultiMappers --outSAMtype $outSAMtype --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart --soloUMIlen $soloUMIlen --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes $outSAMattributes --readFilesIn $R2 $R1";
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
		my $generate_genome = "$EXE --runThreadN $threads --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $fa";
		system("$generate_genome");
	}
}

sub KBref_if_nec {
	
	my $index_KB = "F";
	
	for my $KB_index_file (@KB_index_files) {
		if (-e "$genome_dir/$KB_index_file") {
			if (-s "$genome_dir/$KB_index_file") {
			} else {
				print "$genome_dir/$KB_index_file exists but is empty, will index KB genome\n";
				$index_KB = "T";
				last;
			}
		} else {
			print "Can't find $genome_dir/$KB_index_file, will index KB genome\n";
			$index_KB = "T";
			last;
		}
	}
	
	if ($index_KB eq "F") {
		print "KB genome index files all found, will proceed with the existing index\n";
	} elsif ($index_KB eq "T") {
		my $generate_genome = "$EXE ref -i $genome_dir/transcriptome.idx -g $genome_dir/transcripts_to_genes.txt -f1=$genome_dir/cdna.fa $fa $gtf";
		system("$generate_genome");
	}
}
sub run_KB {
	
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
	
	my $STAR_output_dir = $output_dir . "/alignment_outs/";
	mkdir $STAR_output_dir unless (-d $STAR_output_dir);

	#my $run_STAR = "STAR --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";
	my $run_KBcount = "$EXE count -x=$technology -g=$genome_dir/transcripts_to_genes.txt -i=$genome_dir/transcriptome.idx -o=$output_dir -t=$threads -m=$ram --tmp=~/kbtemp --h5ad --mm $R1 $R2";
	#my $run_KBcount = "$EXE count -x=$technology -g=$genome_dir/transcripts_to_genes.txt -i=$genome_dir/transcriptome.idx -o=$output_dir -t=$threads -m=$ram --h5ad --mm $R1 $R2";


	system("$run_KBcount");
	
	#analyze_BAM($STAR_output_dir);
}
sub convert_h5ad_to_10x {
	my $h5to10x = "$code_dir/python/h5adto10x.py";
	#my $run_STAR = "STAR --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";
	my $run_convert = "python $h5to10x -i $output_dir/counts_unfiltered/adata.h5ad -o $output_dir/alignment_outs/Solo.out/Gene/raw -v $technology";
	system("$run_convert");
	
	#analyze_BAM($STAR_output_dir);
}
sub convert_BUS_to_text {
	#my $run_STAR = "STAR --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";
	my $run_B2t = "bustools text -o $output_dir/output.unfiltered.txt $output_dir/output.unfiltered.bus";
	system("$run_B2t");
	
	#analyze_BAM($STAR_output_dir);
}
sub Alevin_if_nec {
	KBref_if_nec();
	my $cdna_fa = "$genome_dir/cdna.fa";
	my ($genome_name) = $fa =~ /(.*)\.fa/; 
	my $chrnames_txt = "$genome_dir/avchrnames.txt";

	my $gentrome = "$genome_dir/avgentrome.fa";
	my $index = "$genome_dir/avsidx";
	my $tx2gene2 = "$genome_dir/tx2gene.2.txt";

	my @AV_index_files = ($chrnames_txt, $gentrome, $index, $tx2gene2);

	my $index_AV = "F";
	
	for my $AV_index_file (@AV_index_files) {
		if (-e "$AV_index_file") {
			if (-s "$AV_index_file") {
			} else {
				print "$AV_index_file exists but is empty, will index Alevin genome\n";
				$index_AV = "T";
				last;
			}
		} else {
			print "Can't find $AV_index_file, will index Alevin genome\n";
			$index_AV = "T";
			last;
		}
	}

	if ($index_AV eq "F") {
		print "Alevin genome index files all found, will proceed with the existing index\n";
	} elsif ($index_AV eq "T") {
		system("grep \">\" $fa | cut -d \">\" -f 2 | cut -d \" \" -f 1 > $chrnames_txt");
		system("cat $cdna_fa $fa > $gentrome");
		#salmon index
		system("$EXE index -t $gentrome -i $index -p $threads -d $chrnames_txt");
		#transcripts to gene info
		#transcripts_to_genes.txt generated by kb ref
		my $kb_tx2gene = "$genome_dir/transcripts_to_genes.txt";
		system("cut -f1,2 $kb_tx2gene > $tx2gene2");
	}

}
sub run_Alevin {

	my $cdna_fa = "$genome_dir/cdna.fa";
	my ($genome_name) = $fa =~ /(.*)\.fa/; 
	my $chrnames_txt = "$genome_dir/avchrnames.txt";
	my $gentrome = "$genome_dir/avgentrome.fa";
	my $index = "$genome_dir/avsidx";
	my $tx2gene2 = "$genome_dir/tx2gene.2.txt";

	#transcripts_to_genes.txt generated by kb ref
	my $kb_tx2gene = "$genome_dir/transcripts_to_genes.txt";

	for my $R ($R2, $R1) {
		unless ( (-e $R) && (-s $R) ) {
			die "$R is not present or is empty\n";
		}
	}


	#example run
	my $fq1 = $R1;
	my $fq2 = $R2;

	system("$EXE alevin -l ISR -i $index -1 $fq1 -2 $fq2 -o $output_dir -p $threads --tgMap $tx2gene2 --chromium --dumpFeatures --dumpBfh");

	my $final_output_dir = $output_dir . "/alignment_outs/Solo.out/Gene/raw";
	system("mkdir -p $final_output_dir");
	system("mv $output_dir/alevin/quants_tier_mat.gz $final_output_dir/matrix.mtx.gz");
	system("mv $output_dir/alevin/quants_mat_cols.txt $final_output_dir/features.tsv");
	system("mv $output_dir/alevin/quants_mat_rows.txt $final_output_dir/barcodes.tsv");
}
sub CRref_if_nec {
	
	my $index_CR = "F";
	my $CR_ref = $genome_dir . "/cr_references/";

	my $genomeid = "$host_species" . "_host_viruses.$virus_database.with_" . $host_ref_genome;

	if ($virus_database eq "viruSITE.NCBIprokaryotes"){
		$genomeid = "$host_species" . "_host_viruses_microbes.$virus_database.with_" . $host_ref_genome;
	}

	my @CR_index_files = $genomeid; #as of v2.7.5a

	for my $CR_index_file (@CR_index_files) {
		if (-e "$genome_dir/$CR_index_file") {
			if (-s "$genome_dir/$CR_index_file") {
			} else {
				print "$genome_dir/$CR_index_file exists but is empty, will index Cell ranger genome\n";
				$index_CR = "T";
				last;
			}
		} else {
			print "Can't find $genome_dir/$CR_index_file, will index Cell ranger genome\n";
			$index_CR = "T";
			last;
		}
	}
	
	if ($index_CR eq "F") {
		print "Cell ranger genome index files all found, will proceed with the existing index\n";
	} elsif ($index_CR eq "T") {
		my $generate_genome = "$EXE mkref --genome=$genomeid --fasta=$fa --genes=$gtf --nthreads=$threads";

		chdir $genome_dir;
		system("$generate_genome");
		# Move the genome file
		# system("mv $genomeid $genome_dir");

	}

}
sub run_CR {
	my $genomeid = "$host_species" . "_host_viruses.$virus_database.with_" . $host_ref_genome;

	if ($virus_database eq "viruSITE.NCBIprokaryotes"){
		$genomeid = "$host_species" . "_host_viruses_microbes.$virus_database.with_" . $host_ref_genome;
	}
	my $CR_ref = $genome_dir . "/$genomeid";

	for my $R ($R1) {
		unless ( (-e $R) && (-s $R) ) {
			die "$R is not present or is empty\n";
		}
	}
	
	#my $run_STAR = "STAR --runThreadN $threads --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax $max_multi --outFileNamePrefix $STAR_output_dir --sjdbGTFfile $gtf --readFilesCommand $readFilesCommand --soloCBwhitelist $barcodes_whitelist --soloType Droplet --soloBarcodeReadLength $soloBarcodeReadLength --soloStrand $soloStrand --soloUMIfiltering $soloUMIfiltering --soloCellFilter $soloCellFilter --soloCBmatchWLtype 1MM multi pseudocounts --outSAMattributes CR CY CB UR UY UB sM GX GN --readFilesIn $R2 $R1";
	my $run_CR = "$EXE count --id=run_$R2 --fastqs=$R1 --sample=$R2 --localcores=$threads --transcriptome=$CR_ref";

	system("$run_CR");
	#system("mv run_$R2  $genome_dir");

	my $final_output_dir = $output_dir . "/alignment_outs/Solo.out/Gene/raw";

	chdir $output_dir;
	system("mkdir -p $final_output_dir");
	system("mv run_$R2/outs/raw_feature_bc_matrix/* $final_output_dir");
	#system("mv run_$R2 $output_dir");

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
	my $genome_gtf = $genome_dir . "/$host_species" . "_host_viruses.$virus_database.with_" . $host_ref_genome;
	if ($virus_database eq "viruSITE.NCBIprokaryotes"){
		$genome_fa = $genome_dir . "/$host_species" . "_host_viruses_microbes.$virus_database.with_" . $host_ref_genome . ".fa";
		$genome_gtf = $genome_dir . "/$host_species" . "_host_viruses_microbes.$virus_database.with_" . $host_ref_genome;
	}
	unless ( (-e $genome_fa) && (-s $genome_fa) ) {
		die "$genome_fa is not present or is empty\n";
	}
	
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