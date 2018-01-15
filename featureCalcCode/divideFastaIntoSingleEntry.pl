#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use  Bio::Seq;

my $input_file = $ARGV[0] or die "Need to input fasta file on the command line\n";
my $output_dir = $ARGV[1] or die "Need to input output file folder on the command line\n";
my $prefix = $ARGV[2] or die "Need to input prefix on the command line\n";

my $catchseq_seqio_obj = Bio::SeqIO->new(-file=>"$input_file", -format=>'fasta');

my $count = 0;

while(my $seq_obj = $catchseq_seqio_obj->next_seq)
{
    $count++;

    #  在这儿处理每个序列的信息
    my $display_name = $seq_obj->display_name; 
    my $desc = $seq_obj->desc;                                #   序列的描述
	my $seq = $seq_obj->seq;                                   #   序列字符串
	my $seq_type = $seq_obj->alphabet;                    #   序列的类型（dna还是蛋白质？）
	my $seq_length = $seq_obj->length;                    #   序列的长度
    #  以下省略
    
    #print "\$display_name:\n$display_name\n";
    #print "\$desc:\n$desc\n";
    #print "\$seq:\n$seq\n";
    #print "\$seq_type:\n$seq_type\n";
    #print "\$seq_length:\n$seq_length\n";

    my $seqContent = ">".$display_name." ".$desc."\n".$seq."\n";

    print "++++++++++++++++++++++++++++++++++++++++++\n";
    print "Single seq content is as follows: \n";

    print $seqContent;

    print "Writing formated data to $output_dir/$prefix\_$count.fasta \n";
	open (O,">$output_dir/$prefix\_$count.fasta");
	print O $seqContent;
	close O;

	print "success!\n";

	



}

print "solved sequences:\n";
print "$count\n";