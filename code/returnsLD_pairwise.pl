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
#use warnings;

# computes and returns LD values between the given variants (our list, but in a pairwise mode)
#GET ld/:species/pairwise/:id1:/id2/:d_prime/:r2/:population_name

use HTTP::Tiny;
use JSON;
use Data::Dumper;

my $http = HTTP::Tiny->new();
my @arraySNPs = ();

# we open the SNP file of Doritool and store the SNPs in an array
my $VCFfile = $ARGV[0];
my $inputfiletype = ''; # this is to know which type of file it is
if ( $VCFfile =~ /.vcf$/ || $VCFfile =~ /.VCF$/) { $inputfiletype = 'vcf'; }
else { $inputfiletype = 'rs'; }

my $cleanvariant = '';

open (VCFFILE, $VCFfile) or die "Could not open VCF $VCFfile";
while (my $line= <VCFFILE>) {
	chomp $line;
	if ( $line =~ /^#/ || $line =~ /^$/) {next;}
	my @fields = split("\t",$line);
	if ( $inputfiletype eq 'rs') {
	  $cleanvariant = $fields[0];
	  $cleanvariant =~ s/\s*$//g;
	}
	elsif ( $inputfiletype eq 'vcf' ) {
	  $cleanvariant = $fields[2];
	  $cleanvariant =~ s/\s*$//g;
	}
	if ( $cleanvariant !~ /^rs/) {next;}   # this removes the kgp or abc variants (trivial names)
	push(@arraySNPs,$cleanvariant);
}
close(VCFFILE);
print "File type that has been introduced as input: $inputfiletype\n";
print "SNPs for which LD is calculated in a pairwise mode in EUR population (includes CEU, IBS, TSI, FIN, GBR...):\n";
print "@arraySNPs";
print "\n";
;
my $LDfile = "$ENV{'OUTPUTDIR'}/pairwise_LD.tsv";
open (LDFILE, ">$LDfile") or die "Could not open LDFILE $LDfile)";
print LDFILE "pairwise_SNPs\tLD_r2\tLD_Dprime\n";

for (my $snp1=0; $snp1<=$#arraySNPs; $snp1++) {
	for (my $snp2=$snp1+1; $snp2<=$#arraySNPs; $snp2++) {
		my $server = 'http://grch37.rest.ensembl.org';
		my $ext = '/ld/human/pairwise/'.$arraySNPs[$snp1].'/'.$arraySNPs[$snp2].'?population_name=1000GENOMES:phase_3:EUR';
 		print "$arraySNPs[$snp1]\t$arraySNPs[$snp2]\n";
		my $response = $http -> get($server.$ext, { headers => { 'Content-type' => 'application/json' } });
#  		die "Failed!\n" unless $response->{success};   # comment this because stops prematurely the program

		if(length $response->{content}) {
			my $hash = decode_json($response->{content});
			my $claves = keys $hash;   # this is to handle the "error" which is a hash (not array). The error is usually this: 'error' => 'Could not retrieve a variation feature.' And just after this, it launches: hash keys are 1. Not an ARRAY reference at ./returnsLD_pairwise.pl line 67.
#  			print "\nhash keys are $claves\n";
# 			local $Data::Dumper::Terse = 1;
# 			local $Data::Dumper::Indent = 1;
# 			print Dumper $hash;
# 			print "\n";

  			if ( $claves eq 1 ) {
			  if (ref($hash) eq 'HASH') {next;} # forces to leave this because the error are not hash refs or ARRAY, but hashes just as this %hash
  			  else {
# 				print "defined:". ${$hash}[0]-> {'r2'} ;
 				print LDFILE $arraySNPs[$snp1].'/'.$arraySNPs[$snp2]."\t";
				print LDFILE ${$hash}[0] -> {'r2'}."\t".${$hash}[0] -> {'d_prime'}."\n" ;
  			  }
			}
			 else {next;}
		}
	}
}

close(LDFILE);

print "\n";
exit 0;
