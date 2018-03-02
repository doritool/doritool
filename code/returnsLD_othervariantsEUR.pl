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

# computes and returns LD values between the given variant and all other variants in a window centered around the given variant. The window size is set to 500 kb by default and it may be changed.
#GET ld/:species/:id/:d_prime/:r2/:population_name/:window_size

use HTTP::Tiny;
use JSON;
use Data::Dumper;

my $http = HTTP::Tiny->new();
my @arraySNPs = ();

# we open the SNP file of Doritool and store the SNPs in an array
my $VCFfile = $ARGV[0];
my $LDrequired = $ARGV[1];
my $inputfiletype = ''; # this is to know which type of file it is
if ( $VCFfile =~ /.vcf$/ || $VCFfile =~ /.VCF$/) { $inputfiletype = 'vcf'; }
else { $inputfiletype = 'rs'; }

my $cleanvariant = '';

open (VCFFILE, $VCFfile) or die "Could not open VCF $VCFfile";

while (my $line= <VCFFILE>) {
	chomp $line;
	if ( $line =~ /^#/ || $line =~ /^$/) {next;}  # skip header or empty lines
	my @fields = split ("\t",$line);
	  if ( $inputfiletype eq 'rs') {
	  $cleanvariant = $fields[0];
	  $cleanvariant =~ s/\s*$//;
	}
	elsif ( $inputfiletype eq 'vcf' ) {
	  $cleanvariant = $fields[2];
	  $cleanvariant =~ s/\s*$//;
	}
	if ( $cleanvariant !~ /^rs/) {next;}   # this removes the kgp or abc variants (trivial names)
	push(@arraySNPs,$cleanvariant);
}
close(VCFFILE);
print "SNPs for which LD with other variants in a 500 Kb window extracted in EUR population (includes CEU, IBS, TSI, FIN, GBR...): \n";
print "@arraySNPs";
print "\n";

my $LDfile = "$ENV{'OUTPUTDIR'}/other_variants_in_LD_1000GP.tsv";
open (LDFILE, ">$LDfile") or die "Could not open LDFILE $LDfile)";
print LDFILE "variantGWAS\tvariant1000GP\tLD_r2\tLD_Dprime\n";

for (my $snp1=0; $snp1<=$#arraySNPs; $snp1++) {
	my $server = 'http://grch37.rest.ensembl.org';
	my $ext = '/ld/human/'.$arraySNPs[$snp1].'/1000GENOMES:phase_3:EUR?window_size=500;r2='.$LDrequired;
	my $response = $http->get($server.$ext, { headers => { 'Content-type' => 'application/json' }});
# 	die "Failed!\n" unless $response->{success}; # comment this because stops prematurely the program
# 	print $arraySNPs[$snp1]."\n";
	if(length $response->{content}) {
	  my $hash = decode_json($response->{content});
	  local $Data::Dumper::Terse = 1;
	  local $Data::Dumper::Indent = 1;
#  	  print Dumper $hash;
#  	  print "\n";
        my $maxnumber1000G_variants = 10000;
#    	  print "hash type are: $hash\n"; # this is to handle the "error" which is a hash (not array)
 	  if ( ref($hash) eq 'HASH' ) {next;}   # forces to leave this because the error are not hash refs
 	  elsif ( ref($hash) eq 'ARRAY' ) {
		for (my $i=0; $i <= $maxnumber1000G_variants; $i++) {
		    if (defined ${$hash}[$i]) {
			my $variation1 = '';
			my $variation2 = '';
				if (${$hash}[$i] -> {'variation1'} eq $arraySNPs[$snp1]) {
					$variation1 = ${$hash}[$i] -> {'variation1'};
					$variation2 = ${$hash}[$i] -> {'variation2'};
				}
				else {
					$variation1 = ${$hash}[$i] -> {'variation2'};
					$variation2 = ${$hash}[$i] -> {'variation1'};
				}
			print LDFILE $variation1."\t".$variation2."\t".${$hash}[$i] -> {'r2'}."\t".${$hash}[$i] -> {'d_prime'}."\n" ;
		    }
		}
	  }
	  else { next;}
   }
}
print "\n";
exit 0;
