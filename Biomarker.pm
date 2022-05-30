#!/usr/bin/perl
##############################################################################
#
# Biomarker.pm - precision oncology trial and therapy recommender
# 
# Biomarker integration functions
#
# This file is part of Oncology Treatment and Trials Recommender (OTTR) and 
# Precision OTTR (POTTR).
#
# Copyright 2019-2020, Frank Lin & Kinghorn Centre for Clinical Genomics, 
#                      Garvan Institute of Medical Research, Sydney
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
##############################################################################

package Biomarker;

use strict;
use warnings;

use POSIX;

use lib '.';
use TSV;

use Rules;

our @EXPORT;
our @EXPORT_OK;

##################################################################################
BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 1.1 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw( );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}

sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }

my $f_initiailised ;
our %gene_feature_table;
our %biomarker_synonyms;  
our %variant_synonyms;  
our %biomarker_name_lookup = (
	'TMB' => 'tumour_mutational_burden',
	'MMR' => 'mismatch_repair',
	'MSI' => 'microsatellite_instability', 
	'LOH' => 'loss_of_heterozygosity', 
	'HRD' => 'homologous_recombination_deficiency_score', 
);

my %aa_code = (
	'Ala' => 'A',
	'Arg' => 'R', 
	'Asn' => 'N', 
	'Asp' => 'D',
	'Cys' => 'C', 
	'Gln' => 'Q',
	'Glu' => 'E',
	'Gly' => 'G',
	'His' => 'H',
	'Ile' => 'I', 
	'Leu' => 'L', 
	'Lys' => 'K', 
	'Met' => 'M', 
	'Phe' => 'F', 
	'Pro' => 'P', 
	'Ser' => 'S', 
	'Thr' => 'T', 
	'Trp' => 'W', 
	'Tyr' => 'Y', 
	'Val' => 'V',
	'Xaa' => 'X',
	'Unk' => 'X',
	'Asx' => 'B',
	'Glx' => 'Z',
	'Xle' => 'J',
	'Ter' => '*',
);

my $aa_regex = "(?:".join('|', sort keys %aa_code)."|[A-Z])";

$aa_code{$_} = $_ for values %aa_code;


#######################################################################
sub interp_variants {
	&ON_DEMAND_INIT;
	my $biomarker_name = shift;
	my $biomarker_code = shift;
	my $biomarker_spec = shift;

	my @alterations ;
	
	if ( $biomarker_name =~ /\s/ ) {
		$biomarker_name =~ s/\s+/_/g;
		$biomarker_name = lc($biomarker_name);
	} elsif ( exists $biomarker_synonyms{ uc($biomarker_name) } and ( $biomarker_synonyms{uc($biomarker_name)} ne uc($biomarker_name) ) ) {
		$biomarker_name = $biomarker_synonyms{ uc($biomarker_name) } ;
	}
	
	if ( $biomarker_spec =~ /^(?:high|low|proficient|deficient)$/i ) {
		$biomarker_spec = lc($biomarker_spec);
	}
	
	if ( exists $variant_synonyms{ "$biomarker_name:$biomarker_spec" } ) { ; # if synonyms exist
# 		print STDERR "SYNONYM:\t"."$biomarker_name:$biomarker_spec"."\t".$variant_synonyms{ "$biomarker_name:$biomarker_spec" }."\n";
		push @alterations, ( $variant_synonyms{ "$biomarker_name:$biomarker_spec" } =~ s/^$biomarker_name://r ) 
	}
	
	for ( $biomarker_code ) {
		/[VSADEF]/ and do { 
# 			$$facttable{'biomarker'} = $biomarker_name ; 
			push @alterations, 'alteration' ;
		};

		/S/ and do { 
			push @alterations, $biomarker_spec; 
			push @alterations, "splice_mutation"; 
			if ( ($biomarker_spec =~ /exon[[:space:][:punct:]]+(\d+)(?:[[:space:][:punct:]]+(?:skipping|splice))+[[:space:][:punct:]]+(?:variant|mutation)/ ) ) {
				push @alterations, "exon_$1"."_skipping_mutation"; 
				push @alterations, "exon_$1"."_splice_mutation"; 
				push @alterations, "exon_$1"."_deletion"; 
				push @alterations, "exon_$1"."_mutation"; 
			}
			last; 
		};
		
		/E/ and do { push @alterations, $biomarker_spec; last; };
		/A/ and do { push @alterations, ('amplification'); last; };
		/D/ and do { push @alterations, ('deletion', 'homozygous_deletion', 'loss-of-function_mutation'); last; }; # 'truncating_mutation'
		/F/ and do { 
			push @alterations, ('fusion'); 
			if ( ($biomarker_spec =~ /[LR]:([A-Z0-9]+?)--?([A-Z0-9]+)/ ) or #
				 ($biomarker_spec =~ /([A-Z0-9]+?)--?([A-Z0-9]+)[_\s]*(?i:fusions?)/ ) ) {
# 				push @alterations, "$1-$2 fusion" ;
				push @alterations, "$1-$2"."_fusion" ;  # alternative normalised representation with an underscore
			};
			last; 
		};
		
		/T/ and do { # Tumour-wide complex biomarker;
			for ( $biomarker_name ) {
				/(TMB|MMR|MSI|LOH|HRD)/ and do { # FIXME - to replace with a separate dictionary.
					$biomarker_name = $biomarker_name_lookup{$1} if exists $biomarker_name_lookup{$1} ;
					$biomarker_spec =~ s/$1-//;
					$biomarker_spec = lc($biomarker_spec) if $biomarker_spec =~ /High|Low/i;
					last;
				};
			}
			push @alterations, $biomarker_spec;
			last;
		};
		
		/V/ and do {
			my ($aa_org, $aa_pos, $aa_org2, $aa_pos2, $aa_mut, $fs_ter, $ssite_spec, $delins_spec);
			my ($na_org, $na_pos, $na_org2, $na_pos2, $na_alt, $na_mut);
			my @conseq_t; # consequence types
			
			my $f_germline = 0;
			
			for ( $biomarker_spec ) {
				s/(?:,[_ ]?|\()?germline(?:\))?//i and do { $f_germline = 1; };
				
				( ($aa_org, $aa_pos, $aa_mut, $fs_ter) = /(?:p\.)\(??($aa_regex)(\d+)($aa_regex)?fs(\*\d+)?\)?/ ) and do {
					$aa_org = $aa_code{$aa_org};
					$aa_mut = $aa_code{$aa_mut};
					push @conseq_t, 'mutation';
					push @conseq_t, 'frameshift_variant';
					last;
				};
				
				( ($aa_org, $aa_pos, $aa_mut) = /^(?:p\.)?\(?($aa_regex)(\d+)($aa_regex)\)?\)?$/ ) and do {
					$aa_org = $aa_code{$aa_org};
					$aa_mut = $aa_code{$aa_mut};
					push @conseq_t, 'mutation';
					push @conseq_t, 'missense_variant'; # missense_variant
					last;
				};
				
				( ($aa_org, $aa_pos) = /^(?:p\.)?\(?($aa_regex)(\d+)\)?$/ ) and do {
					$aa_org = $aa_code{$aa_org};
					push @conseq_t, 'mutation';
					push @conseq_t, 'missense_variant'; # missense_variant
					last;
				};
				
				( ($aa_org, $aa_pos) = /^(?:p\.)?\(?($aa_regex)(\d+)\*\)?$/ ) and do {
					$aa_org = $aa_code{$aa_org};
					push @conseq_t, 'mutation';
					push @conseq_t, 'stop_gained';
					$aa_mut = '*';
					last;
				};
				
				/splice site/ and ( ($na_org, $ssite_spec) = /(?:c\.)?(\d+)(\-\d+[AGTC]+>[AGTC]+)/ ) and do {
					push @conseq_t, 'mutation';
					push @conseq_t, 'splice_site_variant'; # duplication
					last;
				};
					
				( ($na_pos, $na_pos2, $na_org, $na_alt, $na_mut) = /^(?:c\.)(\d+[\+\-]?\d*)_?(\d+[\+\-]?\d*)?([AGTC]+)?(>|del|delins|ins|dup)([AGTC]+)?/ ) and do {
					push @conseq_t, 'mutation';
					push @conseq_t, 'insertion' if ($na_alt eq 'ins') or ($na_alt eq 'dup') ;
					push @conseq_t, 'indel'  if (defined($na_org) and length($na_org)) != 1 or ( defined($na_mut) and length($na_mut) != 1);
					last;
				};
				
				my $muttype ;
				( ($aa_org, $aa_pos, $aa_org2, $aa_pos2, $muttype, $delins_spec) = /^(?:p\.)?\(?($aa_regex)(\d+)_?($aa_regex)?(\d+)?((?:ins|dup|del(?:ins)?|>))(.*)\)?$/ ) and do {
					$aa_org = $aa_code{$aa_org};
					$aa_org2 = $aa_code{$aa_org2};
					$aa_mut = $aa_code{$aa_mut};
					$muttype = 'delins' if $muttype eq '>';
					push @conseq_t, 'mutation';
					push @conseq_t, 'inframe_insertion' if ($muttype eq 'ins') or ($muttype eq 'dup') ;
					push @conseq_t, 'inframe_deletion'  if ($muttype eq 'del') ;
					push @conseq_t, ('stop_gained', 'truncating_mutation')  if ( $delins_spec eq '*' );
					push @conseq_t, 'indel'  if defined($aa_pos2) and length($delins_spec) != 1 ;
					if ($muttype eq 'delins') {
						my $len_org = ( defined($aa_pos2) ? ($aa_pos2 - $aa_pos + 1) : 1 ) ;
						my $len_mut = length $delins_spec;
						push @conseq_t, 'inframe_indel' ;
						if ( $len_org  < $len_mut ) {
							push @conseq_t, 'inframe_insertion' ;
						} elsif ( $len_org  > $len_mut ) {
							push @conseq_t, 'inframe_deletion' ;
						} 
					}
					last;
				};
			}
			
			push @alterations, $biomarker_spec;
			push @alterations, $aa_org.$aa_pos.$aa_mut  if defined $aa_org and defined $aa_pos and defined $aa_mut ;
			if ( defined $aa_org and defined $aa_pos ) {
				$aa_org = $aa_code{$aa_org};
				for my $conseq_t (@conseq_t) {
					push @alterations, "codon_".$aa_pos."_".$conseq_t; # missense_variant"  ; # and defined $aa_mut 
				}
			}
			
			if ( defined $na_pos and defined $na_alt ) {
				for my $conseq_t (@conseq_t) {
					push @alterations, "nucleotide_".$na_pos."_".$conseq_t; 
				}
			}
			
			for my $conseq_t (@conseq_t) {
				push @alterations, $conseq_t ; 
				push @alterations, 'truncating_mutation' if $conseq_t =~ /frameshift_variant|stop_gained/ ;
			}

			if ( defined $aa_pos ) {
				$aa_pos2 = $aa_pos if ! defined $aa_pos2;
				my @positional_features ;
				for ($biomarker_name) { # Annotate gene features
					exists $gene_feature_table{$biomarker_name} and do {
						for my $row ( @{ $gene_feature_table{$biomarker_name} } ) {
							my $feat_start = $$row{'start'} ;
							my $feat_end   = $$row{'end'} ;
							next if $feat_start =~ /^[cg]\./;
							next if $feat_end =~ /^[cg]\./;
							next if $aa_pos2 < $feat_start ;
							next if $aa_pos  > $feat_end ;
							my $feat = $$row{'feature'};
							push @positional_features, "${feat}_$_" for @conseq_t;
							push @positional_features, "${feat}_mutation";
						}
					};
				}
				
				for my $s (@positional_features) {
					$s =~ s/_variant/_mutation/;
					$s =~ s/_indel/_missense_mutation/;
					$s =~ s/inframe_//;
					push @alterations, $s ;
				}
			}
			
			if ( defined $na_pos ) {
				$na_pos2 = $na_pos if ! defined $na_pos2;
				my ($na_pos_loc, $na_pos_delta)   = ( $na_pos  =~ /(\d+)([\+\-]\d+)?/ );
				my ($na_pos2_loc, $na_pos2_delta) = ( $na_pos2 =~ /(\d+)([\+\-]\d+)?/ );
				
				my @positional_features ;
				for ($biomarker_name) { # Annotate gene features
					exists $gene_feature_table{$biomarker_name} and do {
						for my $row ( @{ $gene_feature_table{$biomarker_name} } ) {
							my $feat_start = $$row{'start'} ;
							my $feat_end   = $$row{'end'} ;
							my $feat       = $$row{'feature'} ;
							
							if ( $feat_start =~ /^[c]\.(\d+)([\+\-]\d+)?/ ) {
								my ($feat_start_loc, $feat_start_delta) = ($1, $2);
								my ($feat_end_loc, $feat_end_delta) = ( $feat_end =~ /^[c]\.(\d+)([\+\-]\d+)?/ );
# 								print ">> $biomarker_name\t$feat_start\t$feat_end\t$na_pos_loc\t$na_pos_delta\t$na_pos2_loc\t$na_pos2_delta\n";
								next if ($na_pos2_loc + ($na_pos2_delta // 0 ) ) < ( $feat_start_loc + ( $feat_start_delta // 0) ) ;
								next if ($na_pos_loc  + ($na_pos_delta  // 0 ) ) > ( $feat_end_loc   + ( $feat_end_delta   // 0) ) ;
							} else { # default is p. 
								next if ( $na_pos2_loc / 3 ) < $feat_start ;
								next if ( $na_pos_loc  / 3 ) > $feat_end ;
							}

							push @positional_features, "${feat}_$_" for @conseq_t;
							push @positional_features, "${feat}_mutation";
						}
					};
				}
				
				for my $s (@positional_features) {
					$s =~ s/_variant/_mutation/;
					$s =~ s/_indel/_missense_mutation/;
					$s =~ s/inframe_//;
					push @alterations, $s ;
				}
			}

			if ( $f_germline ) {
				my @germline_alterations = map { "$_,germline" } @alterations ;
				push @alterations, @germline_alterations ;
			}
			
			last;
		};
	}
	
	my @alterations_uniq ;
	
	for my $a ( @alterations ) {
		next if ! length $a;
		push @alterations_uniq , "$biomarker_name:$a" if ! grep { $_ eq "$biomarker_name:$a" } @alterations_uniq ;
	}
	
	return @alterations_uniq ;
}

#######################################################################
sub load_gene_feature_table {
	for my $srcfile ( POTTRConfig::get_paths('data', 'variant-feature-file') ) {
		for ( file($srcfile) ) { 
			chomp ;
			my ($gene, $feature, $start, $end, $reference) = split /\t/, $_;
			my %a = ( 'start' => $start, 'end' => $end, 'feature' => $feature );
			push @{ $gene_feature_table{$gene} }, \%a;
# 			$gene_feature_table{$gene}{$feature}{'start'} = $start;
# 			$gene_feature_table{$gene}{$feature}{'end'} = $end;
		}
	}
# 	for my $srcfile ( POTTRConfig::get_paths('data', 'biomarker-synonym-file') ) {
# 	}
}

#######################################################################
sub load_biomarker_synonym_file {
	for my $srcfile ( POTTRConfig::get_paths('data', 'biomarker-synonym-file') ) {
		for ( file($srcfile) ) { 
			chomp ;
			my ($a, @b) = split /\t/, $_; 
			for my $m ($a, @b) { 
				$biomarker_synonyms{ uc($m) } = $a;
			}
		}
	}
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'variant-synonym-file') ) {
		for ( file($srcfile) ) { 
			chomp ;
# 			my ($a, @b) = split /\t/, $_; 
# 			for my $m ($a, @b) { 
# 				$variant_synonyms{ $m } = $a;
# 			}
			my (@variants) = split /\t/, $_; 
			for my $x (@variants) { 
				for my $y (@variants) { 
					next if $x eq $y;
					$variant_synonyms{ $x } = $y;
				}
			}
		}
	}
}

#######################################################################
sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;
	load_gene_feature_table;
	load_biomarker_synonym_file ;
}

1;
