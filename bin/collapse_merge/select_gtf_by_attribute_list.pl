#!/usr/bin/env perl
#Author: ALin
#Purpose: To select lines in a gtf file by looking at a list of attribute content in specific attribute.
#Usage: perl select_gtf_by_attribute_list_v1.2.pl <in file> <output> <attribute content> <attribute>
#Change log:
#	v1.2	2021-11	A bug was fixed for handling lines without the specified attribute.

use strict;

if(@ARGV != 4){
	print "#Usage: perl select_gtf_by_attribute_list.pl <in file> <output> <attribute content> <attribute>\n";
	exit;
}

my $in = shift @ARGV;
my $out = shift @ARGV;
my $list = shift @ARGV;
my $attribute = shift @ARGV;

open(IN, "<", $in) or die"Cannot open $in!\n";
open(LIST, "<", $list) or die"Cannot open $list!\n";
open(OUT, ">", $out) or die"Cannot create $out!\n";

my %key_words = ();

while(<LIST>){
	chomp;
	my $key_word = $_;
	$key_word =~ s/\r//g;
	$key_words{$key_word} = 1;
}

close LIST;

my $num = (keys %key_words);
print "Total entries in the list: $num\n";

my $new_line_flag = 0;

while(<IN>){
	my $line = $_;
	chomp $line;
	my @line = split('\t', $line);
	if($line[8] =~ /(^|\s)$attribute \"([^\"]+)\"/){
		my $temp = $2;
		if(exists $key_words{$temp}){
			if($new_line_flag == 0){
				$new_line_flag = 1;
			}
			else{
				print OUT "\n";
			}
			print OUT "$line";
		}
	}
}

close IN;
close OUT;


















