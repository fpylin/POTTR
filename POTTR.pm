#!/usr/bin/perl
##############################################################################
#
# POTTR.pm - precision oncology trial and therapy recommender
# 
# Main inference module
#
# This file is part of Oncology Treatment and Trials Recommender (OTTR) and 
# Precision OTTR (POTTR).
# 
#
# Copyright 2019-2022, Frank Lin & Kinghorn Centre for Clinical Genomics, 
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

use strict;
use warnings;


package POTTR;

use POSIX;

##################################################################################
BEGIN {
	use Exporter   ();
	our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
	$VERSION     = sprintf "%d.%03d", q$Revision: 0.95 $ =~ /(\d+)/g;
	@ISA         = qw(Exporter);
	@EXPORT      = qw();
	@EXPORT_OK   = qw(&score_tier);
}


use lib '.';
use TSV;
use CancerTypes;       
use Rules;
use Biomarker;
use Therapy;
use ClinicalTrials;
use Evidence;
use Data::Dumper;
# use Carp qw(cluck longmess shortmess);


###############################################################################
# Utility functions
###############################################################################

sub min { return undef if (! scalar(@_) ); my $v = shift; for (@_) { $v = $_ if ($_ < $v) ; } return $v; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }
sub trim { my $s = shift; $s =~ s/^\s*|\s*$//g if (defined $s) ; return $s; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub v_safe { my $x = shift; return '' if ! defined $x; return $x; }



##################################################################################
sub score_tier {
	my $x = shift;
	my $sref = shift;
	
	if ( defined($sref) ) { # sref indicates the rank order of the tiers
# 		print "[$x]\t";
		my @a = @{ $sref } ;
		for my $i (0 .. $#a) {
# 			print "$i\n" if lc($a[$i]) eq lc($x);
			return $i if lc($a[$i]) eq lc($x);
		}
		return scalar(@a);
	}
	
	# default using a system similar to OncoKB/TOPOGRAPH
	my ($n) = ( $x =~ /([0-9]+)/ );
	$n  = 5    if ( $x =~ /S/   );
	$n  = 6    if ( $x =~ /U/   );
	$n  = +10  if ( $x =~ /^\(/ );
	$n  = -1   if ( $x =~ /^\(?R1/ );
	$n  = 2.2  if ( $x =~ /^\(?R2/ );
	$n += 0.5  if ( $x =~ /B/   );
	$n += 0.4  if ( $x =~ /R$/  );
	$n  = 5.2  if ( $x =~ /R2B/   );
	$n += 1;
	return $n;
}

##################################################################################
sub is_highest_tier {
	my $x = shift;
	my $sref = shift;
	return 1 if ! defined($sref) and uc($x) eq 'R1' ; # default using a system similar to OncoKB/TOPOGRAPH - R1 is the highest
	return 1 if defined($sref) and lc($$sref[0]) eq lc($x) ; # sref indicates the rank order of the tiers
	return 0;
}


sub max_tier {
	my @tiers = @_;
	my ($max_tier) = sort { score_tier($a) <=> score_tier($b) } grep { defined } @tiers ;
	return $max_tier;
}


#######################################################################
sub new {
    my ($class, $params) = @_;
    my $self = {};
    
    $self->{'params'} = $params ;
    $self->{'modules'} = Rules::Modules->new($params);
    
    my %reason_cache;
    $self->{'reason_cache'} = \%reason_cache;
    
    
    bless($self, $class);
    
    $self->init();
    
    return $self;
}

#######################################################################
# sub get_tier_list {
# 	my $self = shift;
# 	return keys %{ $self->{'tier_list'} }
# }
# 


#######################################################################
sub load_module_init {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('00 - Initialisation and setup');

	$rs->define_dyn_rule( 'define_tier_rank_order',  sub { my $f = shift; my $facts = $_[0];
		my @retval;
		
		for ($f) {
# 			warn "$f";
			/initial-fact/ and do {
				$facts->deffunc('is_resistance_tier', sub { my $facts = $_[0]; my $tier = $_[1];
						if ( my $a = $facts->getvar('tier_list_resistance') ) {
							return ( exists $$a{$tier} ); 		
						} 
						return ( $tier =~ /^R/ );
					}
				);
				last;
			};
			/^tier-rank-order:\s*(.*?)\s*$/ and do {
				my @tier_order = split /\s+/, $1;
				$facts->setvar('tier_order', \@tier_order );
# 				$self->{'tier_order'} = \@tier_order ;
				@retval = ( Facts::mk_fact_str($f, 'processed') );
				last;
			}; 
			/^tiers-list:\s*(.+?)\s*$/ and do {
				my @tier_list = split /\s+/, $1;
				my %tier_list = map { $_ => 1 } @tier_list ;
				$facts->setvar('tier_list', \%tier_list);
				$facts->setvar('tier_list_regex', join("|", ( sort { length($b) <=> length($a) } @tier_list ) ) );
				last;
			};
			/^tiers-list-resistance:\s*(.+?)\s*$/ and do {
				my @tier_list = split /\s+/, $1;
				my %tier_list = map { $_ => 1 } @tier_list ;
				$facts->setvar('tier_list_resistance', \%tier_list);
				last;
			};
		}
		return @retval ;
		}
	);
	
# 	$rs = $self->{'modules'}->add_module('00.2 - Initialisation and setup: predefined rules');

	if ( scalar @POTTRConfig::predefined_rules ) {
# 		print STDERR ( map { "Predefined: $_\n" } @POTTRConfig::predefined_rules ) ;
		$rs->load( @POTTRConfig::predefined_rules ) ;
	}
}


#######################################################################
sub gen_rule_ret_msg { # ($\@)
	my $self = shift;
	my $filename = $_[0];
	my $retval_ref = $_[1];
	$self->{'debug_output'} .= sprintf("$filename: %d rules loaded.\n", scalar(@{ $retval_ref }) );  # printf STDERR 
	return @{ $retval_ref } ;
}


#######################################################################
sub load_module_cancer_type_mapping {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('01A - Cancer type mapping');
	
	$rs->define_dyn_rule( 'translate_catype',  sub { $_=shift; 
			return () if ! /^catype:(.+)/ ;
			my $match = CancerTypes::match_catype_whole_word($1);
			my @retval = ( "(catype-processed)" );
			if ( defined $match ) {
				my @ancerstors = CancerTypes::get_ancestors( $match );
				push @retval, "catype:$_" for (@ancerstors);
				return @retval ;
			}

		return @retval ;
		}
	);
	
# 	$rs = $self->{'modules'}->add_module('01C - Define solid tumours');
	$rs->load(
		mkrule( ["(catype-processed)", "NOT catype:Leukaemia", "NOT catype:Lymphoma"], ["catype:Solid tumour"] ),
		mkrule( ["(catype-processed)", "catype:Leukaemia"], ["catype:Haematologic cancer"] )
	);
}


#######################################################################
sub load_module_variant_feature_mapping {
	my $self = shift;

	my $rs = $self->{'modules'}->add_module('01B - Predefined biomarker rules');
	for my $srcfile ( POTTRConfig::get_paths('data', 'biomarker-rules-file') ) {
		$rs->load( file( $srcfile ) );
	}
	
	$rs = $self->{'modules'}->add_module('01C - Variant Feature Mapping');
	$rs->define_dyn_rule( 'preprocessing_variant',  sub { my $f = shift; my $facts = $_[0];
		my @retval ;
		my @tags = $facts->get_tags($f);
		return () if grep { /\Q(initial-fact)\E/ } @tags ;
		
		# FIXME - To move to a dedicated biomarker ontology database. Translating code to short code
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(TMB|MMR|MSI|LOH|HRD):((?!T:.*))/"$1:".lc($2)/e;
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(?i:TMB|tumour_mutational_burden):(?:T:)?(.*)/"tumour_mutational_burden:".lc($1)/e;
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(?i:MMR|mismatch_repair):(?:T:)?(.*)/"mismatch_repair:".lc($1)/e;
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(?i:MSI|microsatellite_instability):(?:T:)?(.*)/"microsatellite_instability:".lc($1)/e;
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(?i:LOH|loss-of-heterozygosity_score):(?:T:)?(.*)/"loss-of-heterozygosity_score:".lc($1)/e;
		push @retval, Facts::mk_fact_str($f, @tags) if $f =~ s/^(?i:HRD|homologous_recombination_deficiency_score):(?:)?(.*)/"homologous_recombination_deficiency_score:".lc($1)/e;
		
		my ($entity, $etype, $espec) ;
		if ( $f =~ /^(.+?):([VSFTE]):(.+)?$/ ) {
			($entity, $etype, $espec) = ($1, $2, $3);
		} elsif ( $f =~ /^(.+?):([ADF])$/ ) {
			($entity, $etype, $espec) = ($1, $2, '');
		} else {
			for ($f) {
				return @retval if /^(?:catype|tier-rank|has_biomarker|sensitive|resistant|prior_|preferential_|info|\(initial-fact\))/; # FIXME - need a positive list of biomarkers
				/^([^:]+):(.*(?:splice|skipping).*)$/i   and do { ($entity, $etype, $espec) = ($1, 'S', $2); last; };
				/^([^:]+):(?:.*(v[IV]+)(?:\b|$).*)$/i      and do { ($entity, $etype, $espec) = ($1, 'S', $2); last; };
				/^([^:]+):(?:amplification|amplified)$/i       and do { ($entity, $etype, $espec) = ($1, 'A', ''); last; };
				/^([^:]+):(?:deletion|homozygous_deletion)$/i  and do { ($entity, $etype, $espec) = ($1, 'D', ''); last; };
				/^([^:]+):(.*fusion.*)$/i   and do { ($entity, $etype, $espec) = ($1, 'F', $2); last; };
				/^([^:]+):(.*express.*)$/i  and do { ($entity, $etype, $espec) = ($1, 'E', $2); last; };
				/^([^:]+):(.+)$/i           and do { ($entity, $etype, $espec) = ($1, 'V', $2); last; };
				return @retval;
			}
		}
		
		my @alterations = Biomarker::interp_variants( v_safe($entity), v_safe($etype), v_safe($espec) ); 
		my $entity_line = "has_biomarker:".v_safe($entity);
# 		push @tags, "INFERRED:oncogenicity";
		push @retval, ( map { Facts::mk_fact_str($_, @tags) } ($entity_line, @alterations) );
		
		return @retval ;
		}
	);
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'oncogenicity-rules-file') ) {
		$rs->load( file( $srcfile ) );
	}
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'biomarker-rules-file') ) { # biomarker rule file is processed again.
		$rs->load( file( $srcfile ) );
	}
	
	$rs->define_dyn_rule( 'infer_mutation_calls',  sub { my $f = shift; my $facts = $_[0];
		return () if $f !~ /^(\S+?):(?:alteration)(,[_\s]*germline)?\s*$/ ;
		my $germline_suffix = $2 // '';
		my $bm = $1;
		my @tags = grep { /^\Q$bm\E:/ } ( $facts->get_tags($f) );
		return () if ! grep { ! /^INFERRED:|^\Q$bm\E:(?i:oncogenic_mutation|alteration|.*fusion|amplification|.*expression)(?:,[_\s]*germline)?/ } @tags;
# 		return () if $f !~ /^(\S+?):(?:oncogenic_mutation)\s*$/ ;
# 		return ( Facts::mk_fact_str($f, ($facts->get_tags($f)), 'INFERRED:oncogenicity') );
		my $f_inferred_oncogenicity = "$bm:oncogenic_mutation".$germline_suffix ;
		return ( Facts::mk_fact_str( $f_inferred_oncogenicity, 'INFERRED:oncogenicity') ); # ($facts->get_tags($f)), 
		}
	);
}


#######################################################################
sub load_module_clinical_data_module {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('01D - Clinical data module');
	
	$rs->load( $self->gen_rules_drug_db_prior_therapy() ); 
}

#######################################################################
sub load_module_variant_evidence_grading {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('02A - Variant Evidence Grading');

	my $cnt = 0;
	for my $evidence_base_file ( POTTRConfig::get_paths('data', 'therapy-evidence-database-file') ) {
		my ( $evidence_base_label )  = (POTTRConfig::get('therapy-evidence-database-name'))[$cnt] ;
		$evidence_base_label //= "KB$cnt";
		$evidence_base_label =~ s|.*/||;
		$evidence_base_label =~ s|\W.*||;
		$rs->load( Evidence::gen_rule_knowledge_base( $evidence_base_label, $evidence_base_file, $evidence_base_file.".rulescache.txt" ) );
		++ $cnt ;
	}

	$rs = $self->{'modules'}->add_module('02B - Evidence Grade Summary');

	$rs->define_dyn_rule( 'summary_grade',  sub { my $f = shift; my $facts = $_[0]; # my $cache = $_[1];
		return () if $f !~ /^(\S+):(treatment(?:_class)?):(.+?)\s*\(/ ;
		my ($biomarker, $type, $therapy) = ($1, $2, $3);
		my $btt = join(":", $biomarker, $type, $therapy);
# 		return () if exists $$cache{'summary_grade'}{'visited'}{$btt};
# 		$$cache{'summary_grade'}{'visited'}{$btt} = 1;
# 		$cntf{$btt}++; print STDERR "$cache\t$cntf{$btt}\t$btt\n";
		
		my $tier_list_regex = $facts->getvar('tier_list_regex') // '[[:alnum:]]+';

		my @matched_facts ;

# 		if ( exists $$cache{'summary_grade'}{'matched_facts'} ) { # Computationally redundant and expensive search. Using cache
# 			@matched_facts = @{ $$cache{'summary_grade'}{'matched_facts'} } ;
# 		} else {
			@matched_facts = grep { /:$type:.*\(\w+ LOE: ((?:$tier_list_regex))(?:\)|, (?i:inferred|referred) from|, histotype agnostic)/ } ( $facts->get_facts_list() ); # R?[1234S][ABR]?
# 			$$cache{'summary_grade'}{'matched_facts'} = \@matched_facts;
# 		}
		
		my $tier_order_ref = $facts->getvar('tier_order');
		
		# Defining a heuristic to handle variant interactions: resistance and variant interactions.
		# 
		#  Given the scenarios:
		#              Oncogene, S  Oncogene, R  TSG, S (I)
		# Oncogene, R  R            -            
		# TSG, S (I)   S            R            -
		# TSG, R (NI)  S            R            S
		# 
		# TSG = tumour supressor gene
		# A heuristic function for final tiering is:
		#   RS = (TSG1 + TSG2 + TSG3 + ... ) * ( Oncogene1 * Oncogene2 * Oncogene3 )
		# 
		# where 0 = resistant, 1 = sensitising;

		my $has_R_tier = 0;
		
		my @max_tier_tags ;
		
		my %tier_retval = ( max_tier => undef, max_tier_score => undef, ref_biomarker => undef );
		my %tier_TSG = ( max_tier => undef, max_tier_score => undef, ref_biomarker => undef );
		my %tier_TSG_by_gene = ();
		my %tier_treatment_class = ( max_tier => undef, max_tier_score => undef, ref_biomarker => undef );
		my %tier_treatment_class_by_gene = ();
# 		my %tier_drug_class = ( max_tier => undef, max_tier_score => undef, ref_biomarker => undef );

		sub comp_max_tier_and_setvar(\%$$$$) {
			my $var = $_[0];
			my $t   = $_[1];
			my $sc  = $_[2];
			my $bm  = $_[3];
			my $facts = $_[4];
			if ( ! defined $var->{'max_tier'} ) { 
				$var->{'max_tier'} = $t;
				$var->{'max_tier_score'} = $sc;
				$var->{'ref_biomarker_tier'}{$bm} = $t;
			} else { 
				if ( $sc < $var->{'max_tier_score'} ) { 
					$var->{'max_tier'} = $t;
					$var->{'max_tier_score'} = $sc;
				}
				$var->{'ref_biomarker_tier'}{$bm} = max_tier( $var->{'ref_biomarker_tier'}{$bm}, $t );
			}
			
# 			print $facts."\n".longmess( "Callstack:" );
			my $tier_order_ref = $facts->getvar('tier_order');
			


			# ref_biomarker is the output variable
			my @ref_biomarkers = uniq ( ( keys %{ $var->{'ref_biomarker_tier'} } ) );
			
			@ref_biomarkers  = sort { ( score_tier($var->{'ref_biomarker_tier'}{$a}, $tier_order_ref) <=> score_tier($var->{'ref_biomarker_tier'}{$b}, $tier_order_ref) ) || ( $a cmp $b ) }  @ref_biomarkers  ;
			
			my @ref_biomarkers1 = grep { ! ( $facts->func('is_resistance_tier', $var->{'ref_biomarker_tier'}{$_} ) xor $facts->func('is_resistance_tier', $var->{'max_tier'} ) ) } @ref_biomarkers; 
			@ref_biomarkers = @ref_biomarkers1  if scalar @ref_biomarkers1  ;
			
			$var->{'ref_biomarker'} = join (", ", @ref_biomarkers );
# 			print join("\t", $var, $t, $sc, $bm, $var->{'max_tier'}, $var->{'max_tier_score'}, $var->{'ref_biomarker'})."\n";
# 			print Dumper($var)."\n";
		}
		
		sub max_max_tier(\%\%$) {
			my $var1 = $_[0];
			my $var2 = $_[1];
			my $facts = $_[2];
			
			if ( ! defined $var1->{'max_tier'} ) {
				%{$var1} = %{$var2}; 
			}
# 			print STDERR "var1 ".join(" ", map {$_ => $$var1{$_}} grep { $$var1{$_} } sort keys %{$var1}) ."\n";
# 			print STDERR "var2 ".join(" ", map {$_ => $$var2{$_}} grep { $$var2{$_} } sort keys %{$var2}) ."\n";
			
			comp_max_tier_and_setvar(%{$var1}, $var2->{'max_tier'}, $var2->{'max_tier_score'}, $var2->{'ref_biomarker'}, $facts) if defined $var2->{'max_tier'};
		}
		
# 		++ $cntf; print STDERR "! $cntf ! $f\n";		
		for my $fact (@matched_facts) {
			if ( $fact =~ /^(.+?):$type:\Q$therapy\E\s*\(\w+ LOE: ((?:$tier_list_regex))(?:\)|,)/ ) { # R?[1234S][ABUR]?
				my $biomarker = $1; # Get any biomarkers matched to the therapy / therapy class
				my $tier = $2;
				my $is_resistance_tier = $facts->func('is_resistance_tier', $tier);
				$has_R_tier = 1 if $is_resistance_tier ; 
				@max_tier_tags = (@max_tier_tags, ( $facts->get_tags($fact) ));
				
				my $score = score_tier($tier, $tier_order_ref) ;
				if ($type eq 'treatment_class') {
# 					# Treat the "treatment_class" as non-dominant 
					if ( is_highest_tier($tier, $tier_order_ref ) ) { 
						comp_max_tier_and_setvar( %tier_treatment_class, $tier, $score, $biomarker, $facts) ;
					}
					comp_max_tier_and_setvar(%{ $tier_treatment_class_by_gene{$biomarker} }, $tier, $score, $biomarker, $facts) ;
				} elsif ( exists $Evidence::cancer_gene_type{$biomarker} and $Evidence::cancer_gene_type{$biomarker} eq 'TSG' ) {  
					# Reverse resistance/sensitivity ranking when encountering a true tumour suppressor gene, arranged by each gene
# 					print STDERR  "A R $type, $biomarker, $therapy, $tier, $fact\n";
					comp_max_tier_and_setvar(%{ $tier_TSG_by_gene{$biomarker} }, $tier, $score, $biomarker, $facts) ;
# 					comp_max_tier_and_setvar(%{ $tier_TSG_by_gene{$biomarker} }, $tier, $score + ( ( ($type =~ /class/ ) && $is_resistance_tier) ? 1000 : 0 ), $biomarker) ;
				} else { 
# 					print STDERR  "A S $type, $biomarker, $therapy, $tier, $fact\n";
					comp_max_tier_and_setvar(%tier_retval, $tier, $score, $biomarker, $facts) ;
					# If there are a drug class with different sensitivity/resistance profile, sensitivity > resistance;
# 					comp_max_tier_and_setvar(%tier_retval, $tier, $score + ( ( ($type =~ /class/ ) && $is_resistance_tier) ? 1000 : 0 ) , $biomarker) ; 
				}
			}
		}
		

		for my $bm (keys %tier_TSG_by_gene) {
			my $is_resistance_tier = $facts->func('is_resistance_tier', $tier_TSG_by_gene{$bm}{max_tier} );
			my $mt = $tier_TSG_by_gene{$bm}{max_tier};
			my $score = score_tier( $mt, $tier_order_ref) + ( $is_resistance_tier ? 1000 : 0 );
			comp_max_tier_and_setvar(%tier_TSG, $mt, $score, $bm, $facts) ;
		}
		
		for my $bm (keys %tier_treatment_class_by_gene) {
			my $is_resistance_tier = $facts->func('is_resistance_tier', $tier_treatment_class_by_gene{$bm}{max_tier} );
			my $mt = $tier_treatment_class_by_gene{$bm}{max_tier};
			my $score = score_tier( $mt, $tier_order_ref) + ( $is_resistance_tier ? 1000 : 0 );
			comp_max_tier_and_setvar(%tier_treatment_class, $mt, $score, $bm, $facts) ;
		}
		
		if ($type eq 'treatment_class') {
			max_max_tier(%tier_retval, %tier_treatment_class, $facts) ;
		} else {
			max_max_tier(%tier_retval, %tier_TSG, $facts) ;
		}
		
		
		#########################
		my @retval ;
		
		my ($max_tier, $ref_biomarker) = ( $tier_retval{max_tier}, $tier_retval{ref_biomarker} );
		
		$max_tier .= "R" if $has_R_tier and $max_tier !~ /R/; # minor resistance class
		
		if ( defined $max_tier ) {
			my $stem = "recommendation_tier";
			$stem .= '_drug_class' if  $type =~ /class/;
			my $pref_score = score_tier($max_tier, $tier_order_ref);
			my %a = map { $_ => 1 } ( split /\s*,\s*/, $ref_biomarker );
			$ref_biomarker = join(', ', sort keys %a );
			# Added "__untrack__" to the tag to clean up historical reviews
			push @retval, Facts::mk_fact_str("$stem:$therapy:$max_tier", "__untrack__" , "referred_from:$ref_biomarker", "pref_score:$pref_score", @max_tier_tags) ;
		}
		
		return @retval;
		}
	);


	$rs = $self->{'modules'}->add_module('02C - Backpropagating AMP LOE');
	$rs->define_dyn_rule( 'backprop_amp_loe',  sub { my $f = shift; my $facts = $_[0]; # my $cache = $_[1];
		my ($bm, $alt) = split /:/, $f, 2;
		return () if ! defined $alt or index($alt, ':') != -1;
		my @tags = $facts->get_tags($f);
		my @tags_amp_loe = grep { /AMP.LOE/ } @tags ;
		return if ! scalar @tags_amp_loe ;
		my @tags_amp_non_loe = grep { ! /AMP.LOE/ and ( index($_, $bm.":") == 0 ) } @tags ;
# 		print ">> $f [".join(" ", @tags_amp_loe, '|', @tags_amp_non_loe )."]\n";
		for my $ancestor (@tags_amp_non_loe) {
			$facts->assert_tags($ancestor, @tags_amp_loe) if $facts->has_fact($ancestor);
		}
		return ();
		}
	);

	$rs = $self->{'modules'}->add_module('02D - Consolidate AMP LOE');
	$rs->define_dyn_rule( 'consolidate_amp_loe',  sub { my $f = shift; my $facts = $_[0]; # my $cache = $_[1];
		my @tags = $facts->get_tags($f);
		my @tags_amp_loe = grep { /AMP.LOE/ } @tags ;
		return if scalar @tags_amp_loe <= 1;
		@tags_amp_loe = sort @tags_amp_loe ;
		shift @tags_amp_loe ;
		$facts->untag($f, @tags_amp_loe) ;
		return ();
		}
	);

}


#######################################################################
sub load_module_drug_sensitivity_prediction {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('03 - Drug sensitivity prediction from recommendation tiers');

	$rs->define_dyn_rule( 'drug_sensitivity_prediction',  sub { my $f=shift; my $facts = $_[0];
		my $tier_list_regex = $self->{'tier_list_regex'} // '[[:alnum:]]+';
		return () if $f !~ /recommendation_tier(_drug_class)?:(.+?):($tier_list_regex)/;
		my ($type, $therapy, $tier) = ( ($1 // '') , $2, $3); 

		my @tags = $facts->get_tags($f);
		return ( Facts::mk_fact_str("resistant_to$type:$therapy", @tags) ) if $facts->func('is_resistance_tier', $tier) or ( $facts->has_fact("resistant_to$type:$therapy") );
		return ( Facts::mk_fact_str("sensitive_to$type:$therapy", @tags) );
	} );
}


########################################################################
sub gen_rules_drug_db_prior_therapy {
	my $self = shift;
	my @drug_db;
	
	for ( Therapy::get_all_drug_class_pairs() ) {
		my ($d, $dc) = split /\t/, $_;
		my @all_parents = Therapy::get_all_parents($d);
		push @all_parents, "systemic_therapy";
		push @drug_db, mkrule( ["prior_therapy:$d"],  [ ( map { Facts::mk_fact_str("prior_therapy:$_") } @all_parents ) ] );
	}

	for ( Therapy::get_all_drug_class_hierarchy() ) {
		my ($a, $b) = split /\t/, $_;
		push @drug_db, mkrule( ["prior_therapy:$b"], [ Facts::mk_fact_str("prior_therapy:$a") ] );
	}

	return $self->gen_rule_ret_msg( "Prior Therapy:", \@drug_db );
}

#######################################################################
sub gen_rules_drug_db {
	my $self = shift;
	my @drug_db;
	
	for ( Therapy::get_all_drug_class_pairs() ) {
		my ($d, $dc) = split /\t/, $_;
		my $processed_lock = "therapy_sens_predicted:$d";
		
		push @drug_db, mkrule(
			["resistant_to_drug_class:$dc", "NOT $processed_lock", "resistant_to:$d"], 
			[ $processed_lock, Facts::mk_fact_str("resistant_to:$d", "CERTAIN:treatment_drug") ], ['prio:-3']
		);
		push @drug_db, mkrule(
			["resistant_to_drug_class:$dc", "NOT $processed_lock", "NOT resistant_to:$d"], 
			[ $processed_lock, Facts::mk_fact_str("resistant_to:$d", "INFERRED:treatment_drug") ], ['prio:-2']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$dc", "NOT $processed_lock", "NOT resistant_to_drug_class:$dc", "NOT resistant_to:$d", "sensitive_to:$d"], 
			[ $processed_lock, Facts::mk_fact_str("sensitive_to:$d", "CERTAIN:treatment_drug", "CERTAIN:treatment_drug:$d" ) ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$dc", "NOT $processed_lock", "NOT resistant_to_drug_class:$dc", "NOT resistant_to:$d", "NOT sensitive_to:$d"],  
			[ $processed_lock, Facts::mk_fact_str("sensitive_to:$d", "INFERRED:treatment_drug", "INFERRED:treatment_drug:$d") ], ['prio:0']
		);
	}

	for ( Therapy::get_all_drug_class_hierarchy() ) {
		my ($a, $b) = split /\t/, $_;
		push @drug_db, mkrule( 
			["resistant_to_drug_class:$a", "NOT sensitive_to_drug_class:$b"], [ Facts::mk_fact_str("resistant_to_drug_class:$b", "CERTAIN:treatment_drug_class") ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$a", "NOT resistant_to_drug_class:$b"], [ Facts::mk_fact_str("sensitive_to_drug_class:$b", "CERTAIN:treatment_drug_class", "CERTAIN:treatment_drug_class:$b") ]
		);
	}

	for my $combo ( Therapy::get_all_known_combinations() ) {
	
		my $combo_dc_orig = Therapy::get_treatment_class( $combo );
		
# 		print "C\t$combo, $combo_dc_orig\n";
		for my $combo_dc ( Therapy::get_all_matched_offspring_treatment_classes( $combo_dc_orig ) ) {
			next if $combo_dc eq $combo_dc_orig;
			my $processed_lock = "therapy_sens_predicted:$combo_dc";
# 			print "C1\t$combo, $combo_dc\n";
			push @drug_db, mkrule(
				["resistant_to_drug_class:$combo_dc_orig", "NOT $processed_lock", "resistant_to:$combo"], 
				[ $processed_lock, Facts::mk_fact_str("resistant_to_drug_class:$combo_dc", "CERTAIN:treatment_drug_class") ], ['prio:-3']
			);
			push @drug_db, mkrule(
				["resistant_to_drug_class:$combo_dc_orig", "NOT $processed_lock", "NOT resistant_to:$combo"], 
				[ $processed_lock, Facts::mk_fact_str("resistant_to_drug_class:$combo_dc", "INFERRED:treatment_drug_class") ], ['prio:-2']
			);
			push @drug_db, mkrule( 
				["sensitive_to_drug_class:$combo_dc_orig", "NOT $processed_lock", "NOT resistant_to_drug_class:$combo_dc"],  
				[ $processed_lock, Facts::mk_fact_str("sensitive_to_drug_class:$combo_dc", "INFERRED:treatment_drug_class", "INFERRED:treatment_drug_class:$combo_dc") ], ['prio:0']
			);
		}
	}
	
	for my $combo ( Therapy::get_all_known_combinations() ) {
	
		my $combo_dc = Therapy::get_treatment_class( $combo );
		
		my $processed_lock = "therapy_sens_predicted:$combo";
# 		print "D\t$combo, $combo_dc\n";
		push @drug_db, mkrule(
			["resistant_to_drug_class:$combo_dc", "NOT $processed_lock", "resistant_to:$combo"], 
			[ $processed_lock, Facts::mk_fact_str("resistant_to:$combo", "CERTAIN:treatment_drug") ], ['prio:-3']
		);
		push @drug_db, mkrule(
			["resistant_to_drug_class:$combo_dc", "NOT $processed_lock", "NOT resistant_to:$combo"], 
			[ $processed_lock, Facts::mk_fact_str("resistant_to:$combo", "INFERRED:treatment_drug") ], ['prio:-2']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$combo_dc", "NOT $processed_lock", "NOT resistant_to_drug_class:$combo_dc", "NOT resistant_to:$combo", "sensitive_to:$combo"],  
			[ $processed_lock, Facts::mk_fact_str("sensitive_to:$combo", "CERTAIN:treatment_drug", "CERTAIN:treatment_drug:$combo") ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$combo_dc", "NOT $processed_lock", "NOT resistant_to_drug_class:$combo_dc", "NOT resistant_to:$combo", "NOT sensitive_to:$combo"],  
			[ $processed_lock, Facts::mk_fact_str("sensitive_to:$combo", "INFERRED:treatment_drug", "INFERRED:treatment_drug:$combo") ], ['prio:0']
		);
	}
	
	
	return $self->gen_rule_ret_msg( "Therapy:", \@drug_db );
}


#######################################################################
sub load_module_therapy_recomendation_and_grading {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('04A - Existing therapy class grading');

	$rs->define_dyn_rule( 'existing_therapy_grading',  sub { my $f=shift; my $facts = $_[0];
		return () if $f !~ /(resistant|sensitive)_to_drug_class:(.+)/;
		my ($sensitivity, $drug_class) = ($1, $2); 
		$drug_class =~ s/\s*\[.+?\]$//;
		
		my $tier_order_ref = $facts->getvar('tier_order');

		
		my @tags = $facts->get_tags($f);
		my @rec_tiers = map { my $a = $_; $a =~ s/.*:\s*//; $a } grep { /^recommendation_tier/ } @tags ; # search both recommendation_tier and recommendation_tier_drug_class
		@rec_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } @rec_tiers ;
		
		my $group_tier = scalar(@rec_tiers) ? $rec_tiers[0] : ( $sensitivity =~ /resistant/i ? 'R' : 'S' );
		for my $t ( grep { /^recommendation_tier:/ } @tags ) { # importing all tags from matched therapies
			next if ! $facts->has_fact($t);
			push @tags, ( $facts->get_tags($t) );
		}
		return ( Facts::mk_fact_str("recommendation_tier_drug_class:$drug_class:$group_tier", @tags) );
		}
	);

	$rs = $self->{'modules'}->add_module('04B - New therapy recommendation and grading (sensitive)');
# 	
	$rs->load( $self->gen_rules_drug_db() );
# 	
}

#######################################################################
sub load_module_therapy_prioritisation { # including all therapies
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('05 - Therapy prioritisation');
	
	$rs->define_dyn_rule( 'therapy_prioritisation',  sub { my $f=shift; my $facts = $_[0];
		my $tier_order_ref = $facts->getvar('tier_order');
		
		return () if $f !~ /(sensitive)_to:(.+)/; # (_drug_class)?
		my ($sens_resi, $therapy) = ($1, $2, $3); # $is_drug_class, 
		my @tags = $facts->get_tags($f);
		my @tiers;
		my @tiers_dc;
		
		for my $tag ( @tags ) { # importing all tags from matched therapies
			next if $tag !~ /^recommendation_tier(_drug_class)?:.+:([[:alnum:]]+)$/ ;
			my ($is_drug_class_rt, $tier) = ($1, $2);
			if ($is_drug_class_rt) {
				push @tiers_dc, $tier;
			} else {
				push @tiers, $tier;
			}
		}
		
		my $max_tier_dc = max_tier(@tiers_dc, 'U');
		my $max_tier    = max_tier(@tiers,    'U');
		
# 		print "\e[1;43m!!!$therapy\t".join("/", $max_tier, $max_tier_dc)."\e[0m\n";
		push @tags, "therapy_recommendation_tier_therapy:$max_tier";
		push @tags, "therapy_recommendation_tier_drug_class:$max_tier_dc";

		my $score = # Score the therapy 
			( pow(20,1) * ( score_tier( $max_tier,    $tier_order_ref ) ) ) + 
			( pow(20,0) * ( score_tier( $max_tier_dc, $tier_order_ref ) ) );
			
		push @tags, "therapy_recommendation_score:$score";
		
		return ( Facts::mk_fact_str("therapy_recommendation:$therapy:$max_tier", @tags) );
	} );
}



#######################################################################
sub load_module_preferential_trial_matching {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('06 - Preferential trial matching');
	
	my @trial_database_files = POTTRConfig::get_paths('data', 'clinical-trial-database-file'); 

	# multiple rule files are allowed to filter search results trials based on clinical eligibility
	my @trial_eligibility_file = POTTRConfig::get_paths('data', 'clinical-trial-eligibility-file'); 
	
	for my $trial_database_file (@trial_database_files) {
		$rs->load( $self->gen_rule_ret_msg(
				$trial_database_file, [
				 ClinicalTrials::gen_rules_clinical_trials(
					$trial_database_file,
					@trial_eligibility_file, 
					$trial_database_file.".rulescache.txt"
				) ]
			)
		);
	}
}

#######################################################################
sub load_module_preferential_trial_prioritisation {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('06 - Preferential trial prioritisation');
	$rs->define_dyn_rule( 'clinical_trial_grading',  sub { my $f=shift; my $facts = $_[0];
		return () if $f !~ /^(?:preferential_trial_id:)(\S+)/;
		my $trial_id = $1; 
		my @tags = $facts->get_tags($f);
# 		my @inf_drugs        = map { my $a = $_; $a =~ s/(?:CERTAIN|INFERRED):treatment_drug://; $a }       grep { /(?:CERTAIN|INFERRED):treatment_drug:/ }       @tags ;
# 		my @inf_drug_classes = map { my $a = $_; $a =~ s/(?:CERTAIN|INFERRED):treatment_drug_class://; $a } grep { /(?:CERTAIN|INFERRED):treatment_drug_class:/ } @tags ;
		my @inf_drugs        = map { /recommendation_tier:(.+):\S+/; $1 }            grep { /recommendation_tier:.+:\S+/ }       @tags ;
		my @inf_drug_classes = map { /recommendation_tier_drug_class:(.+):\S+/; $1 } grep { /recommendation_tier_drug_class:.+:\S+/ } @tags ;
		
		my @trial_match_criteria = map { s/trial_match_criteria://; $_ } grep { /^trial_match_criteria:/} @tags;
		my ($trial_phase)    = map { s/^phase://; $_ } grep { /^phase:/} @tags;
		my @rec_drugs_tier;
		
		my @facts_list = $facts->get_facts_list();
		my @new_tags;
# 		print "TTT $trial_id\t$_\n" for sort @tags;
# 		print "!!! $trial_id\t$_\n" for @inf_drug_classes;
# 		print "!!? $trial_id\t$_\n" for grep { /(?:CERTAIN|INFERRED):treatment_drug_class:/ } @tags ;
# 		push @new_tags, ( grep { /(?:CERTAIN|INFERRED):treatment_drug/ } @tags );
		
		my $tier_order_ref = $facts->getvar('tier_order');

		my %a;
		
		
		
		my @transitive_efficacy_tiers = ('U');       # = maximum tier of the referring recommendation tier.
		
		FACT1: for my $f ( grep { /^recommendation_tier:/} @facts_list) {
			for my $d (@inf_drugs) {
				$f =~ /^recommendation_tier:\Q$d\E:(\S+)/ and do { push @transitive_efficacy_tiers, $1; }  ; # next FACT1 
			}
		}
		
		my @transitive_class_efficacy_tiers = ('U'); # = maximum tier of the referring recommendation tier.
		FACT2: for my $f ( grep { /^recommendation_tier_drug_class:/} @facts_list) {
			for my $dc (@inf_drug_classes) {
				$f =~ /^recommendation_tier_drug_class:\Q$dc\E:(\S+)/ and do { push @transitive_class_efficacy_tiers, $1; }  ; # next FACT2 
			}
		}

		my @biomarker_tiers = ('U');
		FACT3: for my $t ( grep { /^AMP.LOE:/} @tags) {
			$t =~ /^AMP.LOE:(\S+)/ and do { push @biomarker_tiers, $1; next FACT3 }  ;
		}
		
		@transitive_efficacy_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } map { Therapy::tidy_tier($_) } @transitive_efficacy_tiers ;
		@transitive_class_efficacy_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } map { Therapy::tidy_tier($_) } @transitive_class_efficacy_tiers ;
		@biomarker_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } map { Therapy::tidy_tier($_) } @biomarker_tiers;
		
		my $transitive_efficacy_tier = $transitive_efficacy_tiers[0];
		my $transitive_class_efficacy_tier = $transitive_class_efficacy_tiers[0];
		my $biomarker_tier = $biomarker_tiers[0];
		
# 		my $trial_match_criteria_score = (
# 			( ( grep { $_ eq 'drug_class_sensitivity' } @trial_match_criteria ) ? 2 : 0 ) +
# # 			( ( grep { $_ eq 'drug_sensitivity' } @trial_match_criteria ) ? 2 : 0 ) +
# 			( ( grep { $_ eq 'cancer_type' } @trial_match_criteria ) ? 2 : 0 ) 
# 			);
		my $trial_match_criteria_score = (
			( ( ( grep { $_ eq 'drug_class_sensitivity' } @trial_match_criteria ) || ( grep { $_ eq 'drug_sensitivity' } @trial_match_criteria ) ) ? 2 : 0 ) +
# 			( ( grep { $_ eq 'drug_sensitivity' } @trial_match_criteria ) ? 2 : 0 ) +
			( ( grep { $_ eq 'cancer_type' } @trial_match_criteria ) ? 2 : 0 ) 
			);
		my $trial_drug_sensitivity_score = ( ( grep { $_ eq 'drug_sensitivity' } @trial_match_criteria ) ? 2 : 0 ) ;
		
		my $num_referring_drug_classes_score = scalar(@inf_drug_classes) ;
# 		
		push @new_tags, # Annotate the trial
			"transitive_class_efficacy_tier:". ( $a{'transitive_class_efficacy'}  = $transitive_class_efficacy_tier ) ,
			"transitive_efficacy_tier:".       ( $a{'transitive_efficacy'}        = $transitive_efficacy_tier ) ,
			"drug_maturity_tier:".             ( $a{'drug_maturity'}              = ClinicalTrials::get_drug_maturity_tier_by_trial_id( $trial_id, @inf_drugs ) ) , 
			"drug_class_maturity_tier:".       ( $a{'drug_class_maturity'}        = ClinicalTrials::get_drug_class_maturity_tier_by_trial_id( $trial_id, @inf_drug_classes ) ) ,
			"combo_maturity_tier:".            ( $a{'combo_maturity'}             = ClinicalTrials::get_drug_combination_maturity_tier_by_trial_id( $trial_id, @inf_drugs ) ) ,
			"combo_class_maturity_tier:".      ( $a{'combo_class_maturity'}       = ClinicalTrials::get_drug_class_combination_maturity_tier_by_trial_id( $trial_id, @inf_drug_classes ) ),
			"trial_phase_tier:".               ( $a{'trial_phase_tier'}           = ClinicalTrials::get_tier_by_phases_of_trials( $trial_id, $trial_phase ) ),
			"biomarker_tier:".                 ( $a{'biomarker_tier'}             = $biomarker_tier ),
			"trial_match_criteria_score:".     $trial_match_criteria_score ,
			"referring_drug_classes_score:".   $num_referring_drug_classes_score
# 			"transitive_class_efficacy_tiers_str:".join("-", @transitive_class_efficacy_tiers)
		;
		
		
		my $score = # Score the trial
			( pow(20,9) * ( 8 - $trial_match_criteria_score ) ) + 
			( pow(20,8) * ( score_tier( $a{'transitive_class_efficacy'}, $tier_order_ref) ) ) + 
			( pow(20,7) * ( 8 - $trial_drug_sensitivity_score ) ) + 
			( pow(20,6) * ( score_tier( $a{'transitive_efficacy'}, $tier_order_ref) ) ) + 
			( pow(20,5) * ( (19 - $num_referring_drug_classes_score) ) ) +
			( pow(20,4) * ( score_tier( $a{'drug_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,3) * ( score_tier( $a{'trial_phase_tier'}, $tier_order_ref ) ) ) + 
			( pow(20,2) * ( score_tier( $a{'drug_class_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,1) * ( score_tier( $a{'combo_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,0) * ( score_tier( $a{'combo_class_maturity'}, $tier_order_ref ) ) ) +
			0;
			
		push @new_tags, "pref_trial_score:$score";
		
		return ( Facts::mk_fact_str("preferential_trial_id:$trial_id", @new_tags) ) ;
	} );
}

#######################################################################
sub load_module_deinit {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('07A - Tidying up results');
	$rs->define_dyn_rule( 'tidy_inference_tracking',  sub { my $f=shift; my $facts = $_[0];
		my @tags = $facts->get_tags($f);
		my %type;
		for (@tags) {
			/^(CERTAIN|INFERRED):(.*)/ and do { push @{ $type{$2} }, $1; };
		}
		for my $t ( keys %type ) {
			my @s = sort @{ $type{$t} };
			shift @s; # Keep the first item
			$facts->untag($f, "$_:$t") for @s;
		}
		return ();
	} );

	$rs = $self->{'modules'}->add_module('07B - Remove non-dominant contingency assertions');
	$rs->define_dyn_rule( 'remove_contingency_assertions',  sub { my $f=shift; my $facts = $_[0];
		my @tags = $facts->get_tags($f);
		my @completely_matched_items = grep { defined } map { /preferential_trial_complete_match:(.+)/ ? $1 : undef } @tags ;
		return () if ! scalar @completely_matched_items ;
		
		for my $item (@completely_matched_items) {
			my @tags_to_remove = grep { /^\s*\*$item:/ } @tags;
			my %tags_to_remove = map { $_ => 1 } @tags_to_remove ;
			next if ! scalar @tags_to_remove;
			$facts->untag($f, @tags_to_remove);
			# $facts->assert_tags($f, "*Contingency constraint satisified by one path");
		}
		return ();
	} ) ;

	$rs = $self->{'modules'}->add_module('07C - Calculate LOM');
	$rs->define_dyn_rule( 'calculate_LOM',  sub { my $f=shift; my $facts = $_[0];
		my @tags = $facts->get_tags($f);
		my %certain;
#		my %inferred;
		
		for (@tags) {
			/^(?:CERTAIN):([^:]*)/  and do { my $c = $1; $c =~ s/_drug//; $certain{$c} += 1; next ; };
			/^(?:INFERRED):([^:]*)/ and do { my $c = $1; $c =~ s/_drug//; $certain{$c} += 0; next ; }; # $c =~ s/_class$//; 
		}
		
		my %LOM_string = ( 1 => 'Complete match', 2 => 'Partial match', 3 => 'Completely inferred' );
		my $LOM = 0;
		
		if ( scalar(keys %certain) ) {
			$LOM += ( scalar(grep { $certain{$_} == 0 } (keys %certain)) > 0 ) ? 2 : 1;
			$LOM += ( scalar(grep { $certain{$_} > 0 } (keys %certain)) == 0 ) ? 1 : 0;
		}
 		
		if ($LOM) {
			my $LOM_string = "LOM:$LOM - $LOM_string{$LOM}";
			$LOM_string.= join("/", " - LOMreason", (map { "$certain{$_}-$_" } keys %certain) ); #."\n";
			return Facts::mk_fact_str($f, $LOM_string);
		} 
		return ();
	} );

}

#######################################################################
sub init {
	my $self = shift;
	
    my $params = $self->{'params'} ;
    
	$params->{'modules'} = join(',', qw(CTM CLI VFM VEG DSO TRG PTM PTP)) if ! defined $params->{'modules'} ;

	my %modules = map { $_ => 1 } split( /\s*[;, ]\s*/, $params->{'modules'} );
    
# 	print STDERR ">>2 ".$POTTRConfig::f_initialised."\n";
# 	print STDERR ">>2 $_.\n" for @POTTRConfig::predefined_rules;
	POTTRConfig::ON_DEMAND_INIT();
    $self->load_module_init();
	exists $modules{'CTM'}  and $self->load_module_cancer_type_mapping();
	exists $modules{'CLI'}  and $self->load_module_clinical_data_module();
	exists $modules{'VFM'}  and $self->load_module_variant_feature_mapping();
	exists $modules{'VEG'}  and $self->load_module_variant_evidence_grading();
	exists $modules{'DSO'}  and $self->load_module_drug_sensitivity_prediction();
	exists $modules{'TRG'}  and $self->load_module_therapy_recomendation_and_grading();
	exists $modules{'TRG'}  and $self->load_module_therapy_prioritisation();
	exists $modules{'PTM'}  and $self->load_module_preferential_trial_matching();
	exists $modules{'PTP'}  and $self->load_module_preferential_trial_prioritisation();
	$self->load_module_deinit();
}

#######################################################################
sub reason {
	my $self = shift;
	my @facts = @_;
	my $fact_str = join("\n", @facts);
	my @retval ;

	if ( exists $self->{'reason_cache'}{$fact_str} ) {
		@retval = @{ $self->{'reason_cache'}{$fact_str} };
	} else {
		@retval = $self->{'modules'}->run(@facts);
		$self->{'reason_cache'}{$fact_str} = \@retval ;
	}

# 	$self->{'modules'}->print();
	
	$self->{'debug_output'} = $self->{'modules'}{'debug_output'};
	return @retval ;
}


1;
