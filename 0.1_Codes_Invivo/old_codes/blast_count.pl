#!/usr/bin/perl
####
# Assumes output generated by blast with:
# -outfmt "6 qacc sacc evalue qstart qend sstart send"
# Counts total occurrences of top hit for each result.
####

%reads = ();
%counts = ();

# Initialize counts to 0 from reference file.
open (INFILE, "./blastdb/blastdb.fa") or die "Reference file not found.";
while ($row = <INFILE>) {
	if ($row =~ /^>/) {
		$row =~ s/\n//g;
		$row =~ s/\r//g;
		$row =~ s/>//g;
		$counts{$row} = 0;
	}
}
close (INFILE);

open (INFILE, $ARGV[0]) or die "$ARGV[0] not found.";
while ($row = <INFILE>) {
	@rowAry = split (/\t/, $row);
	$rowAry[1] =~ s/\s//g;
	if (not defined $reads{$rowAry[0]}) { # Only count the top (first) hit for each read.
		$reads{$rowAry[0]} = 1;
		if (not defined $counts{$rowAry[1]}) {
			$counts{$rowAry[1]} = 1;
		} else {
			$counts{$rowAry[1]}++;
		}
	}
}
close (INFILE);

for my $key ( sort keys %counts ) {
	print $key . "\t" . $counts{$key} . "\n";
}


