#!/usr/bin/env perl
=start
Author: ALin
Purpose: To parse a gtf file of transcript and prepare exonic regions, or 3' exonic regions.
=cut

use strict;
use Getopt::Long;

my $usage = "Usage: perl prepare_intrapriming_rescue_bed.pl
	-i <String> Input gtf
	-o <String> Output basename
	-d <Integer> The distance from the 3' of the terminal exon (Default: 1000)
	-h <Boolean> Help
";

my ($in, $out, $dist, $help) = ("", "", 1000, 0);

GetOptions(
	'i=s'	=>	\$in,
	'o=s'	=>	\$out,
	'd=i'	=>	\$dist,
        'h!'    =>      \$help,
);

if($help){
	print $usage;
	exit;
}

unless($in || $out){
	print "$usage";
	exit;
}

my $transcript_ref = read_gtf($in);
my ($all_ref, $terminal_ref) = parse_transcript($transcript_ref, $dist);
print_bed($all_ref, $terminal_ref, $out);

sub read_gtf{
	my ($temp_in) = @_;
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_in |") or die "<ERROR> Cannot open $temp_in!\n";
	my %temp_transcript = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[2] ne "exon"){
			next;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $1;
		}
		else{
			print STDERR "<ERROR> transcript_id missing in $temp_in at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
		}
		my @temp_exon = ($temp_line[3] - 1, $temp_line[4]);
		push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
	}
	close $temp_fh;

	return \%temp_transcript;
}

sub parse_transcript{
	my ($temp_transcript_ref, $temp_dist) = @_;
	my %temp_all = ();
	my %temp_terminal = ();
	foreach my $temp_transcript_id (keys %{$temp_transcript_ref}){
		my ($temp_chr, $temp_strand) = ($temp_transcript_ref->{$temp_transcript_id}{'chr'}, $temp_transcript_ref->{$temp_transcript_id}{'strand'});
		unless(exists $temp_all{$temp_chr}){
			%{$temp_all{$temp_chr}} = ();
			%{$temp_terminal{$temp_chr}} = ();
		}
		unless(exists $temp_all{$temp_chr}{$temp_strand}){
			@{$temp_all{$temp_chr}{$temp_strand}} = ();
			@{$temp_terminal{$temp_chr}{$temp_strand}} = ();
		}
		my $temp_num_exon = @{$temp_transcript_ref->{$temp_transcript_id}{'exon'}};
		for(my $i = 0; $i < $temp_num_exon; $i++){
			my @temp_exon = @{$temp_transcript_ref->{$temp_transcript_id}{'exon'}[$i]};
			if($temp_strand eq '-'){
				if($i == 0){
					$temp_exon[0] -= $temp_dist;
					push(@{$temp_terminal{$temp_chr}{$temp_strand}}, \@temp_exon);
				}
			}
			else{
				if($i == ($temp_num_exon - 1)){
					$temp_exon[1] += $temp_dist;
					push(@{$temp_terminal{$temp_chr}{$temp_strand}}, \@temp_exon);
				}
			}
			push(@{$temp_all{$temp_chr}{$temp_strand}}, \@temp_exon);
		}
	}

	return (\%temp_all, \%temp_terminal);
}

sub print_bed{
	my ($temp_all_ref, $temp_terminal_ref, $temp_out) = @_;
	open(my $temp_all_fh, "| sort -k1,1 -k2,2n > ${temp_out}_all.bed") or die "Cannot create ${temp_out}_all.bed!\n";
	open(my $temp_terminal_fh, "| sort -k1,1 -k2,2n > ${temp_out}_terminal.bed") or die "Cannot create ${temp_out}_terminal.bed!\n";
	foreach my $temp_chr (keys %$temp_all_ref){
		foreach my $temp_strand (keys %{$temp_all_ref->{$temp_chr}}){
			my @temp_exon = interval($temp_all_ref->{$temp_chr}{$temp_strand});
			foreach my $temp_exon_ref (@temp_exon){
				print $temp_all_fh "$temp_chr\t$temp_exon_ref->[0]\t$temp_exon_ref->[1]\t${temp_chr}_$temp_exon_ref->[0]_$temp_exon_ref->[1]_${temp_strand}\t0\t$temp_strand\n";
			}
		}
	}
	close $temp_all_fh;
	foreach my $temp_chr (keys %$temp_terminal_ref){
		foreach my $temp_strand (keys %{$temp_terminal_ref->{$temp_chr}}){
			my @temp_exon = interval($temp_terminal_ref->{$temp_chr}{$temp_strand});
			foreach my $temp_exon_ref (@temp_exon){
				print $temp_terminal_fh "$temp_chr\t$temp_exon_ref->[0]\t$temp_exon_ref->[1]\t${temp_chr}_$temp_exon_ref->[0]_$temp_exon_ref->[1]_${temp_strand}\t0\t$temp_strand\n";
			}
		}
	}
	close $temp_terminal_fh;

	return 1;
}

sub interval{
        my ($sref,$start) = @_;
        return if (ref($sref) ne 'ARRAY');

        if (!defined $start) {
                if (wantarray) {
                        my @tmpsets = map {[@{$_}]} @{$sref};
                        $sref = \@tmpsets;
                }
                @{$sref} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$sref};
                $start = 0;
        }
        my $last = $sref->[$start];
        ++$start;

        if (@{$last}){
                for my $ndx ($start .. @{$sref}-1){
                        my $cur = $sref->[$ndx];
                        next if (!@{$cur});

                        if ($cur->[0] >= $last->[0] && $cur->[0] <= $last->[1] ){
                                $last->[1] = $cur->[1] if ($cur->[1] > $last->[1]);
                                @{$cur} = ();
                        }
                        else{
                                last;
                        }
                }
        }
        interval($sref, $start) if ( $start < @{$sref});
        if(wantarray){
                return sort {$a->[0] <=> $b->[0]} map {@{$_} ? $_ : () } @{$sref};
        }
}




