#!/usr/bin/perl
use strict; use warnings;

my $genome_fa = "$ref_dir/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.fa";
#cdna fa is the one generated by kb ref
my $cdna_fa = "$ref_dir/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.cdna.fa";

my ($genome_name) = $genome_fa =~ /(.*)\.fa/; 

my $genome_gtf = "$ref_dir/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.removed_amb_viral_exon.gtf";

my $chrnames_txt = "$genome_name.chrnames.txt";
system("grep \">\" $genome_fa | cut -d \">\" -f 2 | cut -d \" \" -f 1 > $chrnames_txt");

#concatenate cdna fa and genome fa - not mentioned in minnow tutorial but mentioned in Alevin tutorial (https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
my $gentrome = $genome_name . ".gentrome.fa";
system("cat $cdna_fa $genome_fa > $gentrome");

#salmon index
my $index = "$genome_name.sidx";
system("salmon index -t $gentrome -i $index -p $threads -d $chrnames_txt");

#transcripts to gene info
my $tx2gene2 = "$genome_name.tx2gene.2.txt";
#transcripts_to_genes.txt generated by kb ref
my $kb_tx2gene = "$ref_dir/human_host_viruses_microbes.viruSITE.NCBIprokaryotes.with_hg38.transcripts_to_genes.txt";
system("cut -f1,2 $kb_tx2gene > $tx2gene2");

#example run
my $fq1 = "/storage/holab/gastric_cancer_organoid_scRNA-seq/GX001-TO-5Seq_R1.fastq.gz";
my $fq2 = "/storage/holab/gastric_cancer_organoid_scRNA-seq/GX001-TO-5Seq_R2.fastq.gz";

system("salmon alevin -l ISR -i $index -1 $fq1 -2 $fq2 -o $alevin_out_dir -p $threads --tgMap $tx2gene2 --chromium --dumpFeatures --dumpBfh");
