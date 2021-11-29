use strict;
use warnings;

if (@ARGV < 1) {die "Usage: parse_XMLs.pl Metadata.json";}


my $file = $ARGV[0];
open(my $ifh, $file) or die $!;
my $ID = "NONE";
while(<$ifh>) {
	if ($_ =~ /"file_name": "(.+)\.htseq\.counts\.gz"/) {
		print "\n$1";
	}
	if ($_ =~ /"case_id": "(.+)"/) {
		print "\t$1";
	}
	if ($_ =~ /"tumor_stage": (.+)/) {
		print "\tstage:$1";
	}
	if ($_ =~ /"days_to_death": (.+)/) {
		print "\tdeath:$1";
	}
	if ($_ =~ /"days_to_recurrence": (.+)/) {
		print "\trecurr:$1";
	}
	if ($_ =~ /"number_proliferating_cells": (.+)/) {
		print "\tprolif:$1";
	}
	if ($_ =~ /"percent_necrosis": (.+)/) {
		print "\tnecro:$1";
	}
	if ($_ =~ /"percent_normal_cells": (.+)/) {
		print "\tnorm:$1";
	}
	if ($_ =~ /"vital_status": (.+)/) {
		print "\tvital:$1";
	}
	if ($_ =~ /"days_to_last_follow_up": (.+)/) {
		print "\tfollowup:$1";
	}
} close ($ifh);
