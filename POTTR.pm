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

use strict;
use warnings;


package POTTR;

use POSIX;

##################################################################################
BEGIN {
	use Exporter   ();
	our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
	$VERSION     = sprintf "%d.%03d", q$Revision: 0.9 $ =~ /(\d+)/g;
	@ISA         = qw(Exporter);
	@EXPORT      = qw();
	@EXPORT_OK   = qw(&score_tier);
}


use lib '.';
use TSV;
use DOID;       # qw( DO_match_catype DO_name);
use Rules;
use Biomarker;
use Therapy;
use ClinicalTrials;
use Evidence;

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
	
	# If a tier ranking not defined, use OncoKB tiering system.
	my ($n) = ( $x =~ /([0-9]+)/ );
	$n  = 5    if ( $x =~ /S/   );
	$n  = 6    if ( $x =~ /U/   );
	$n  = +10  if ( $x =~ /^\(/ );
	$n  = -1   if ( $x =~ /^\(?R1/ );
	$n  = 2.2  if ( $x =~ /^\(?R2/ );
	$n += 0.4  if ( $x =~ /B/   );
	$n += 0.5  if ( $x =~ /R$/  );
	$n += 1;
	return $n;
}

sub max_tier {
	my @tiers = @_;
	my ($max_tier) = sort { score_tier($a) <=> score_tier($b) } @tiers ;
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
sub get_tier_list {
	my $self = shift;
	return keys %{ $self->{'tier_list'} }
}

sub is_resistance_tier {
	my $self = shift;
	my $tier = shift;
	if ( exists $self->{'tier_list_resistance'} ) {
		return ( exists $self->{'tier_list_resistance'}{$tier} ); 
	} else {
		return ( $tier =~ /^R/ );
	}
}


#######################################################################
sub load_module_init {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('00 - Initialisation and setup');

	$rs->load( @POTTRConfig::predefined_rules ) if scalar @POTTRConfig::predefined_rules ;

	$rs->define_dyn_rule( 'define_tier_rank_order',  sub { my $f = shift; my $facts = $_[0];
		my @retval;
		for ($f) {
			/^tier-rank-order:\s*(.*?)\s*$/ and do {
				my @tier_order = split /\s+/, $1;
				$self->{'tier_order'} = \@tier_order ;
				@retval = ( Facts::mk_fact_str($f, 'processed') );
				last;
			}; 
			/^tiers-list:\s*(.+?)\s*$/ and do {
				my @tier_list = split /\s+/, $1;
				$self->{'tier_list'}{$_} ++ for @tier_list ;
				$self->{'tier_list_regex'} = join("|", ( sort { length($b) <=> length($a) } @tier_list ) );
				last;
			};
			/^tiers-list-resistance:\s*(.+?)\s*$/ and do {
				my @tier_list = split /\s+/, $1;
				$self->{'tier_list_resistance'}{$_} ++ for @tier_list ;
				last;
			};
		}
		return @retval ;
		}
	);
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
		my ($match) = DOID::DO_match_catype($1);
# 		print "C\t$1\t$match\n";
# 		print join("\n", map { "$_\t$DO_name{$_}" } DOID::get_ancestors($match))."\n" if defined $match;
		my @retval = ( "(catype-processed)" );
		if ( defined($match) and exists($DO_name{$match}) ) {
			my @ancerstors = DOID::get_ancestors( $match );
			for my $a (@ancerstors) {
				push @retval, ( "catype:".$a, "catype:".$DO_name{$a}, "catype_name:$DO_name{$a}" ) ;
			}
			return @retval ;
		}

		return @retval ;
		}
	);
	
# 	$rs = $self->{'modules'}->add_module('01C - Define solid tumours');
	$rs->load(
		mkrule( ["(catype-processed)", "NOT catype:leukemia"], ["catype:Solid Tumour"] ),
		mkrule( ["(catype-processed)", "catype:leukemia"], ["catype:Liquid Cancer"] )
	);
}


#######################################################################
sub load_module_variant_feature_mapping {
	my $self = shift;

	my $rs = $self->{'modules'}->add_module('01B - Variant Feature Mapping');
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
				return @retval if /^(?:catype|tier-rank|has_biomarker|sensitive|resistant|prior_|\(initial-fact\))/; # FIXME - need a list of biomarkers
				/^([^:]+):(?:amplification|amplified)$/i       and do { ($entity, $etype, $espec) = ($1, 'A', ''); last; };
				/^([^:]+):(?:deletion|homozygous_deletion)$/i  and do { ($entity, $etype, $espec) = ($1, 'D', ''); last; };
				/^([^:]+):(.*fusion.*)$/i   and do { ($entity, $etype, $espec) = ($1, 'F', $2); last; };
				/^([^:]+):(.*express.*)$/i  and do { ($entity, $etype, $espec) = ($1, 'E', $2); last; };
				/^([^:]+):(.*splice.*)$/i   and do { ($entity, $etype, $espec) = ($1, 'S', $2); last; };
				/^([^:]+):(.+)$/i           and do { ($entity, $etype, $espec) = ($1, 'V', $2); last; };
				return @retval;
			}
		}
		
		my @alterations = Biomarker::interp_variants( v_safe($entity), v_safe($etype), v_safe($espec) ); 
		my $entity_line = "has_biomarker:".v_safe($entity);
		push @retval, ( map { Facts::mk_fact_str($_, @tags) } ($entity_line, @alterations) );
		
		return @retval ;
		}
	);
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'oncogenicity-rules-file') ) {
		$rs->load( file( $srcfile ) );
	}
	
	$rs->define_dyn_rule( 'infer_mutation_calls',  sub { my $f = shift; my $facts = $_[0];
		return () if $f !~ /^(\S+?):(?:alteration)(,[_\s]*germline)?\s*$/ ;
		my $germline_suffix = $2 // '';
# 		return () if $f !~ /^(\S+?):(?:oncogenic_mutation)\s*$/ ;
# 		return ( Facts::mk_fact_str($f, ($facts->get_tags($f)), 'INFERRED:oncogenicity') );
		my $f_inferred_oncogenicity = "$1:oncogenic_mutation".$germline_suffix ;
		return ( Facts::mk_fact_str( $f_inferred_oncogenicity, 'INFERRED:oncogenicity') ); # ($facts->get_tags($f)), 
		}
	);
}


#######################################################################
sub load_module_clinical_data_module {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('01C - Clinical data module');
	
# 	$rs->load( $self->gen_rules_drug_db_prior_therapy() ); # Not currently active.
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
# 		my ($evidence_base_label) = ( $evidence_base_file =~ /^(?:.+\/)?(\w+)/ ); # "GIMRKB";
		$rs->load( Evidence::gen_rule_knowledge_base( $evidence_base_label, $evidence_base_file, $evidence_base_file.".rulescache.txt" ) );
		++ $cnt ;
	}

	
	$rs = $self->{'modules'}->add_module('02B - Evidence Grade Summary');
	$rs->define_dyn_rule( 'summary_grade',  sub { my $f = shift; my $facts = $_[0];
		return () if $f !~ /^(\S+):(treatment(?:_class)?):(.+?)\s*\(/ ;
		my ($biomarker, $type, $therapy) = ($1, $2, $3);
		my $ref_biomarker = undef; # gene referred from
		my $max_tier = undef;
		my $max_tier_score = undef;
		my @max_tier_tags ;
		my $tier_list_regex = $self->{'tier_list_regex'} // '[[:alnum:]]+';

		my @matched_facts = grep { /:$type:.*\(\w+ LOE: ((?:$tier_list_regex))(?:\)|, referred from|, histotype agnostic)/ } ( $facts->get_facts_list() ); # R?[1234S][ABR]?
		
		my $tier_order_ref = $self->{'tier_order'} ;
		
		my $has_R_tier = 0;
		for my $fact (@matched_facts) {
			if ( $fact =~ /^.+?:$type:\Q$therapy\E\s*\(\w+ LOE: ((?:$tier_list_regex))(?:\)|,)/ ) { # R?[1234S][ABUR]?
				my $tier = $1;
				$has_R_tier = 1 if $self->is_resistance_tier($tier); 
# 				print "A $type, $biomarker, $therapy, $tier\n";
				my $score = score_tier($tier, $tier_order_ref) ;
				if ( ( ( ! defined $max_tier ) or ( $score < $max_tier_score ) ) ) { 
					$max_tier = $tier ;
					$max_tier_score = $score ; 
					$ref_biomarker = $biomarker;
					@max_tier_tags = $facts->get_tags($fact);
				} elsif ( $score == $max_tier_score ) { 
					$ref_biomarker .= ", $biomarker" if ! $ref_biomarker =~ /$biomarker/;
					@max_tier_tags = (@max_tier_tags, ( $facts->get_tags($fact) ));
				}
			}
		}
		
		$max_tier .= "R" if $has_R_tier and $max_tier !~ /R/; # minor resistance class
		
		my @retval ;
		
		if ( defined $max_tier ) {
			my $stem = "recommendation_tier";
			$stem .= '_drug_class' if  $type =~ /class/;
			my $pref_score = score_tier($max_tier, $tier_order_ref);
			push @retval, Facts::mk_fact_str("$stem:$therapy:$max_tier", "__untrack__", "referred_from:$ref_biomarker", "pref_score:$pref_score", @max_tier_tags) ;
		}

		return @retval;
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
		return ( Facts::mk_fact_str("resistant_to$type:$therapy", @tags) ) if $self->is_resistance_tier($tier) or ( $facts->has_fact("resistant_to$type:$therapy") );
		return ( Facts::mk_fact_str("sensitive_to$type:$therapy", @tags) );
	} );
}


# #######################################################################
# sub gen_rules_drug_db_prior_therapy {
# 	my $self = shift;
# 	my @drug_db;
# 	
# 	for ( Therapy::get_all_drug_class_pairs() ) {
# 		my ($d, $dc) = split /\t/, $_;
# 		push @drug_db, mkrule( ["prior_therapy:$d"],  [ Facts::mk_fact_str("prior_therapy:$dc") ] );
# 	}
# 
# 	for ( Therapy::get_all_drug_class_hierarchy() ) {
# 		my ($a, $b) = split /\t/, $_;
# 		push @drug_db, mkrule( ["prior_therapy:$b"], [ Facts::mk_fact_str("prior_therapy:$a") ] );
# 	}
# 
# 	return $self->gen_rule_ret_msg( "Prior Therapy:", \@drug_db );
# }

#######################################################################
sub gen_rules_drug_db {
	my $self = shift;
	my @drug_db;
	
	for ( Therapy::get_all_drug_class_pairs() ) {
		my ($d, $dc) = split /\t/, $_;
		push @drug_db, mkrule(
			["resistant_to_drug_class:$dc"], [ Facts::mk_fact_str("resistant_to:$d", "INFERRED:treatment_drug") ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$dc", "NOT resistant_to_drug_class:$dc", "NOT resistant_to:$d"],  [ Facts::mk_fact_str("sensitive_to:$d", "INFERRED:treatment_drug") ]
		);
	}

	for ( Therapy::get_all_drug_class_hierarchy() ) {
		my ($a, $b) = split /\t/, $_;
		push @drug_db, mkrule( 
			["resistant_to_drug_class:$a", "NOT sensitive_to_drug_class:$b"], [ Facts::mk_fact_str("resistant_to_drug_class:$b", "CERTAIN:treatment_drug_class") ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$a", "NOT resistant_to_drug_class:$b"], [ Facts::mk_fact_str("sensitive_to_drug_class:$b", "CERTAIN:treatment_drug_class") ]
		);
	}

	for my $combo ( Therapy::get_all_known_combinations() ) {
		my $combo_dc = Therapy::get_treatment_class( $combo );
# 		print "C $combo, $combo_dc\n";
		push @drug_db, mkrule(
			["resistant_to_drug_class:$combo_dc"], [ Facts::mk_fact_str("resistant_to:$combo", "INFERRED:treatment_drug") ], ['prio:-1']
		);
		push @drug_db, mkrule( 
			["sensitive_to_drug_class:$combo_dc", "NOT resistant_to_drug_class:$combo_dc", "NOT resistant_to:$combo"],  [ Facts::mk_fact_str("sensitive_to:$combo", "INFERRED:treatment_drug") ]
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
		
		my $tier_order_ref = $self->{'tier_order'} ;
		
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
		my $tier_order_ref = $self->{'tier_order'} ;
		
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

	my $trial_eligibility_file = POTTRConfig::get_first_path('data', 'clinical-trial-eligibility-file'); 
	
	for my $trial_database_file (@trial_database_files) {
		$rs->load( $self->gen_rule_ret_msg(
				$trial_database_file, [
				 ClinicalTrials::gen_rules_clinical_trials(
					$trial_database_file,
					$trial_eligibility_file, 
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
		my @annotations;
		my ($trial_id, $tags) = ($1, $2); 
		my @tags = $facts->get_tags($f);
		my @inf_drugs        = map { s/INFERRED:treatment_drug://; $_ }       grep { /INFERRED:treatment_drug:/ }       @tags ;
		my @inf_drug_classes = map { s/INFERRED:treatment_drug_class://; $_ } grep { /INFERRED:treatment_drug_class:/ } @tags ;
		my @rec_drugs_tier;
		
		my @facts_list = $facts->get_facts_list();
		my @new_tags;
		
		my $tier_order_ref = $self->{'tier_order'} ;

		my %a;
		
		my @transitive_efficacy_tiers = ('U');       # = maximum tier of the referring recommendation tier.
		
		
		FACT1: for my $f ( grep { /^recommendation_tier:/} @facts_list) {
			for my $d (@inf_drugs) {
				$f =~ /^recommendation_tier:\Q$d\E:(\S+)/ and do { push @transitive_efficacy_tiers, $1; next FACT1 }  ;
			}
		}
		
		my @transitive_class_efficacy_tiers = ('U'); # = maximum tier of the referring recommendation tier.
		FACT2: for my $f ( grep { /^recommendation_tier_drug_class:/} @facts_list) {
			for my $dc (@inf_drug_classes) {
				$f =~ /^recommendation_tier_drug_class:\Q$dc\E:(\S+)/ and do { push @transitive_class_efficacy_tiers, $1; next FACT2 }  ;
			}
		}
		
		@transitive_efficacy_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } map { Therapy::tidy_tier($_) } @transitive_efficacy_tiers ;
		@transitive_class_efficacy_tiers = sort { score_tier($a, $tier_order_ref) <=> score_tier($b, $tier_order_ref) } map { Therapy::tidy_tier($_) } @transitive_class_efficacy_tiers ;
		
		my $transitive_efficacy_tier = $transitive_efficacy_tiers[0];
		my $transitive_class_efficacy_tier = $transitive_class_efficacy_tiers[0];
		
		push @new_tags, # Annotate the trial
			"transitive_class_efficacy_tier:". ( $a{'transitive_class_efficacy'}  = $transitive_class_efficacy_tier ) ,
			"transitive_efficacy_tier:".       ( $a{'transitive_efficacy'}        = $transitive_efficacy_tier ) ,
			"drug_maturity_tier:".             ( $a{'drug_maturity'}              = ClinicalTrials::get_drug_maturity_tier_by_trial_id( $trial_id, @inf_drugs ) ) , 
			"drug_class_maturity_tier:".       ( $a{'drug_class_maturity'}        = ClinicalTrials::get_drug_class_maturity_tier_by_trial_id( $trial_id, @inf_drug_classes ) ) ,
			"combo_maturity_tier:".            ( $a{'combo_maturity'}             = ClinicalTrials::get_drug_combination_maturity_tier_by_trial_id( $trial_id, @inf_drugs ) ) ,
			"combo_class_maturity_tier:".      ( $a{'combo_class_maturity'}       = ClinicalTrials::get_drug_class_combination_maturity_tier_by_trial_id( $trial_id, @inf_drug_classes ) )
		;
		
		my $score = # Score the trial
			( pow(20,4) * ( score_tier( $a{'transitive_efficacy'}, $tier_order_ref) ) ) + 
			( pow(20,5) * ( score_tier( $a{'transitive_class_efficacy'}, $tier_order_ref) ) ) + 
			( pow(20,3) * ( score_tier( $a{'drug_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,2) * ( score_tier( $a{'drug_class_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,1) * ( score_tier( $a{'combo_maturity'}, $tier_order_ref ) ) ) + 
			( pow(20,0) * ( score_tier( $a{'combo_class_maturity'}, $tier_order_ref ) ) );
			
		push @new_tags, "pref_trial_score:$score";
		
		return ( ( Facts::mk_fact_str("preferential_trial_id:$trial_id", @new_tags) ), @annotations );
	} );
}

#######################################################################
sub load_module_deinit {
	my $self = shift;
	my $rs = $self->{'modules'}->add_module('07 - Tidying up results');
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


	$rs = $self->{'modules'}->add_module('07B - Calculate LOM');
	$rs->define_dyn_rule( 'calculate_LOM',  sub { my $f=shift; my $facts = $_[0];
		my @tags = $facts->get_tags($f);
		my @certain;
		my @inferred;
		
		for (@tags) {
			/^(?:CERTAIN):(.*)/  and do { push @certain,  $1; next ; };
			/^(?:INFERRED):(.*)/ and do { push @inferred, $1; next ; };
		}
		
		my %LOM_string = ( 1 => 'Complete match', 2 => 'Partial match', 3 => 'Completely inferred' );
		my $LOM = 0;
		
		if ( scalar(@inferred) + scalar(@certain) ) {
			$LOM += ( scalar(@inferred) > 0 ) ? 2 : 1;
			$LOM += ( scalar(@certain) == 0 ) ? 1 : 0;
		}
		
		if ($LOM) {
			my $LOM_string = "LOM:$LOM - $LOM_string{$LOM}";
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