use strict;
use warnings;

if (@ARGV < 1) {die "Usage: parse_XMLs.pl clinical.xml";}

my %FIELDS = ("days_to_death"=>1, "days_to_last_followup"=>1, "days_to_last_known_alive"=>1, "clinical_stage"=>1, "pathologic_stage"=>1, "days_to_new_tumor_event_after_initial_treatment"=>1, "vascular_tumor_cell_type"=>1, "patient_id"=>1, "bcr_patient_barcode"=>1, "bcr_patient_uuid"=>1, "file_uuid"=>1);

print(join("\t", keys(%FIELDS))."\n");
foreach my $file (@ARGV) {

	my %data= ();
	open(my $ifh, $file) or die $!;
	my $ID = "NONE";
	while(<$ifh>) {
		if ($_ =~ /<(.+)>(.+)<\/(.+):(.+)>/) {
			my $field = $4;
			my $val = $2;
			if ($field eq "file_uuid"){
				$ID=$val;
			}
			$data{$field} = $val;
		}
	} close ($ifh);
	print $ID;
	foreach my $key (keys(%FIELDS)) {
		if (exists($data{$key})) {
			print "\t$data{$key}";
		} else {
			print "\tNA";
		}
	}
	print "\n";
}
