#!/usr/bin/env perl
=start
Author: ALin
Purpose: To merge two tables by adding table two to table one based on an identifier column.
Usage: perl merge2tables_${version}.pl <table 1> <table 1 identifier column (0-based)> <table 2> <table 2 identifier column> <out> <shift mode> <ID column>
Change log:
	v1.2	2021-03	One argument was added to specify whether to shift the columns or not if no match was found.
	v1.3	2022-02	Added the ID column arguement to indicate whether to add the ID column of the second table.
	v1.4	2022-09	Changed the way to remove the ID column.
	v1.5	2023-01 Added options and enabled case insensitive mode.
	v1.6	2023-05	Fixed a bug when handling supplementary table containing empty cells.
=cut

use strict;
use Getopt::Long;

my $version = "v1.6";

my $usage = "Usage: perl merge2tables.pl
	-t1 <String> Main table
	-c1 <Integer> Main table 1 identifier column (0-based)
	-t2 <String> Supplementary table
	-c2 <Integer> Supplementary table identifier column (0-based)
	-o <String> Output
	-s <Boolean> Shift mode
	-d <Boolean> Ouput the second ID column
	-n <Boolean> Case insensitive
	-h <Boolean> Help
";

my ($t1, $c1, $t2, $c2, $out, $shift, $id_col, $case, $help) = ("", 0, "", 0, "", 0, 0, 0, 0, 0);

GetOptions(
	't1=s'	=>	\$t1,
	'c1=s'	=>	\$c1,
	't2=s'	=>	\$t2,
	'c2=s'	=>	\$c2,
	'o=s'	=>	\$out,
	's!'	=>	\$shift,
	'd!'	=>	\$id_col,
	'n!'	=>	\$case,
	'h!'	=>	\$help
);

unless($t1 && $t2 && $out){
	print STDERR "$usage";
	exit;
}


my ($t2_ref, $num_t2_col) = read_supplementary($t2, $c2, $id_col, $case);
read_main($t2_ref, $num_t2_col, $t1, $c1, $shift, $case, $out);

sub read_supplementary{
	my ($temp_t2, $temp_c2, $temp_id_col, $temp_case) = @_;
	my %temp_t2 = ();
	my $temp_num_t2_col = 0;
	my $temp_fh;
	open($temp_fh, $temp_t2) or die "Cannot open $temp_t2!\n";
	while(<$temp_fh>){
        	chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//g;
		my @temp_line = split('\t', $temp_line, -1);
		if($temp_num_t2_col < @temp_line){
			$temp_num_t2_col = @temp_line;
		}
		my $temp_id;
		if($temp_case){
			$temp_id = uc($temp_line[$temp_c2]);
		}
		else{
			$temp_id = $temp_line[$temp_c2];
		}
		unless($temp_id_col){
			my @temp_new_line = ();
			for(my $i = 0; $i < @temp_line; $i++){
				if($i != $temp_c2){
					push(@temp_new_line, $temp_line[$i]);
				}
			}
			$temp_line = join("\t", @temp_new_line);
		}
		$temp_t2{$temp_id} = $temp_line;
	}
	close $temp_fh;
	unless($temp_id_col){
		$temp_num_t2_col--;
	}

	return (\%temp_t2, $temp_num_t2_col);
}

sub read_main{
	my ($temp_t2_ref, $temp_num_t2_col, $temp_t1, $temp_c1, $temp_shift, $temp_case, $temp_out) = @_;
	my ($temp_in_fh, $temp_out_fh);
	open($temp_in_fh, $temp_t1) or die "Cannot open $temp_t1!\n";
	open($temp_out_fh, "> $temp_out") or die "Cannot create $temp_out\n";
	while(<$temp_in_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//g;
		my @temp_line = split('\t', $temp_line);
		print $temp_out_fh "$temp_line";
		if($temp_case){
			if(exists $temp_t2_ref->{uc($temp_line[$temp_c1])}){
				print $temp_out_fh "\t$temp_t2_ref->{uc($temp_line[$temp_c1])}";
			}
			else{
				if($temp_shift){
					for(my $i = 0; $i < $temp_num_t2_col; $i++){
						print $temp_out_fh "\t";
					}
				}
			}
		}
		else{
			if(exists $temp_t2_ref->{$temp_line[$temp_c1]}){
				print $temp_out_fh "\t$temp_t2_ref->{$temp_line[$temp_c1]}";
			}
			else{
				if($temp_shift){
					for(my $i = 0; $i < $temp_num_t2_col; $i++){
						print $temp_out_fh "\t";
					}
				}
			}
		}
		print $temp_out_fh "\n";
	}
	close $temp_in_fh;
	close $temp_out_fh;

	return 1;	
}







