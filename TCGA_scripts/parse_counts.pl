use strict;
use warnings;

if (@ARGV < 1) {die "Usage: parse_counts.pl .htseq.counts";}

my %ExprMat = ();
my @IDs = ();

foreach my $file (@ARGV) {
	$file =~ /(.+)\.htseq\.counts/;
	my $ID = $1;
	push(@IDs, $ID);
	open(my $ifh, $file) or die $!;
	while(<$ifh>) {
		my @stuff = split(/\s+/);
		$stuff[0] =~ /(ENSG\d+)/;
		my $gene = $1;
		$ExprMat{$gene}->{$ID} = $stuff[1];
	} close ($ifh);
}

print(join("\t", @IDs)."\n");
foreach my $g (sort(keys(%ExprMat))){
	print "$g";
	foreach my $id (@IDs) {
		if (exists($ExprMat{$g}->{$id})) {
			print "\t".$ExprMat{$g}->{$id};
		} else {
			print "\tNA";
		}
	}
	print "\n";
}

