#!/usr/bin/env perl

# DoriTool is the most complete bioinformatic tool offering functional 'in silico' annotation of variants
# Copyright (C) 2017 Isabel Martín-Antoniano, Lola Alonso,Miguel Madrid; Evangelina López de Maturana; Núria Malats. Genetic and Molecular Epidemiology Group of Spanish National Cancer Research Centre (CNIO)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.



use strict;
use warnings;
#use Data::Dumper;

# Take out from start and end regions (./output/VEP/condeloutputclean0.txt) in order to search cytobands (cytoband.tsv) and the variants which overlap those regions as well. output file: condeloutputclean1.txt
# Usage: ./cytoband_finder.pl ./output/VEP/condeloutputclean0.txt ./databases/cytoBand.tsv
# Usage: ./cytoband_finder.pl "$WORKDIR"/VEP/condeloutputclean0.txt ../databases/cytoBand.tsv
my $positionsfile = $ARGV[0];
my $cytobandfile = $ARGV[1];

my $summaryfile = $positionsfile;
$summaryfile =~ s/condeloutputclean0.txt$/condeloutputclean1.txt/;

#1- Read cytoBand file (ENSEMBL Format but without chrprefix). Create a hash output with the following fields: chr, start, end, cytoband. Upload in a hash. The 5th field (GIEMSA positive or negative, and intensity) are not interesting. Beware with the assembly patches like HSCHR4, etc
my %cytobands = ();	# is a hash table
open (CYTOFILE, $cytobandfile) or die "Could not open CYTO file $cytobandfile";
while (my $line = <CYTOFILE>) {
	my @fields = split ("\t",$line);
	$cytobands{$fields[0]}{$fields[1].'_'.$fields[2]} = $fields[3]; # cytobands{chr}{start_end} = cytoband
}
close (CYTOFILE);

# print Dumper keys %cytobands;
# print Dumper keys $cytobands{'16'};
# print Dumper $cytobands{'16'}{'28100000_34600000'};

#3- Read condeloutputclean0.txt file. Take out the position and insert in the last field the cytoband column.
open (TABLEFILE, $positionsfile) or die "Could not open TABLE file $positionsfile";
open (SUMFILE, ">$summaryfile") or die "Could not open SUM file $summaryfile";
while (my $line = <TABLEFILE>) {
	chomp $line;
	if ($line =~ /^Uploaded/) { print SUMFILE "$line\tcytoband\n"; next;}
	my @fields = split ("\t",$line);
	my $variantposition = $fields[1];
	my $chr = $variantposition;
	$chr =~ s/:.*$// ;  # this removes the coordinates also for the HSCHR_ patches and MCHs
	my $coordinate = $variantposition;
	$coordinate =~ s/^.*:// ;
	my $cytoband = '';
	if ( $chr !~ /^[0-9]/ ) {
	  $cytoband = '';  # because it will paste it to chr in line 55
	}
	else {
	  foreach my $startend (keys $cytobands{$chr}) {
		no warnings 'numeric';
		my $cytostart = $startend;  $cytostart =~ s/_.*$//; 	# print "$cytostart\t";
		my $cytoend = $startend;  	$cytoend =~ s/^.*_//;	# print "$cytoend\n";
		if ( $coordinate >= $cytostart && $coordinate <= $cytoend ) {
		   $cytoband .= $cytobands{$chr}{$startend};
		}
		else {next;}
	  }
	}
	$cytoband = $chr.$cytoband;
	print SUMFILE "$line\t$cytoband\n";
}
close (TABLEFILE);
close (SUMFILE);
exit;
