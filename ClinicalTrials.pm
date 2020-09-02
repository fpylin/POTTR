#!/usr/bin/perl
##############################################################################
#
# ClinicalTrials.pm - precision oncology trial and therapy recommender
# 
# Clinical trials assessment functions
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


package ClinicalTrials;

use strict;
use warnings;

use POSIX;

use lib '.';
use DOID;
use Rules;
use Therapy;
use TSV;

our @EXPORT;
our @EXPORT_OK;

##################################################################################
BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 1.1 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}

sub mtime { my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($_[0]); return $mtime; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }


##################################################################################

our $f_initiailised = undef;

my %drug_trial_phases_drugs ; # indexed;
my %drug_trial_phases_drug_classes ; # indexed;
my %drug_trial_phases_combo ; # indexed;
my %drug_trial_phases_combo_classes ; # indexed;


my %drug_pubmed_count ; # = map { my ($d, $c, @rest) = split /\t/, $_; Therapy::get_preferred_drug_name($d) => $c } file( 'data/drug_pubmed_count.txt' );


sub get_trial_href($) {
	my $x = shift;
	for ($x) {
		m|^(NCT\d+)$| and do { return "https://clinicaltrials.gov/ct2/show/$1" };
		m|^(ACTRN\d+)$| and do { return "https://www.anzctr.org.au/TrialSearch.aspx#&&searchTxt=$1" } ;
	}
	return "javascript:void(0)";
}


####################################################################################################
sub write_array_uniq($\@) {
	my $outfile = $_[0];
	my @data = @{ $_[1] };
	my %visited;
	open OUTF, ">$outfile" or die "FATAL: Unable to write $outfile: $!";
	for (@data) {
		next if exists $visited{$_} ;
		$visited{$_} = 1;
		print OUTF "$_\n";
	}
	close OUTF;
	chmod 0666, $outfile;
# 	print $!;
}

####################################################################################################
sub get_tier_by_phases_of_trials {
	&ON_DEMAND_INIT ;
	my @phase_list = @_;
# 	print join(" ", @phase_list)."\n";
	my $max_phase = 0;
	for my $ph ( @phase_list ) {
		next if $ph =~ /Not Applicable/i;
		next if $ph !~ /Phase/i;
		$ph =~ s/(?:Phase|Early)\s*//g;
		$ph =~ s/.*\///g;
		$max_phase = ( $ph gt $max_phase ) ? $ph : $max_phase ;
	}
	
	$max_phase =~ tr/01234/U4321/;
	my $tier = Therapy::tidy_tier($max_phase);
	return $tier;
}

####################################################################################################
my %drug_maturity_tier_cache;
sub get_drug_maturity_tier {
	&ON_DEMAND_INIT ;
	my $trial_id ;
	$trial_id = shift if ( $_[0] =~ /NCT|ACTRN/ );
	my @inf_drugs = @_;
	my @phases;
	
	for my $d (@inf_drugs) {
		$d = Therapy::get_preferred_drug_name($d) ;
		if ( exists $drug_maturity_tier_cache{$d} ) {
			push @phases, $drug_maturity_tier_cache{$d};
			next;
		}
		my @phases_drug ;
		if ( Therapy::is_a( $d, 'regulatory_approved' ) ) {
			push @phases_drug, 'Phase 4';
			next;
		} elsif ( exists( $drug_trial_phases_drugs{$d} ) ) {
			for my $trial ( @{ $drug_trial_phases_drugs{$d} } ) { # @drug_trial_phases_drugs
				next if defined $trial_id and ($$trial{'trial_id'} eq $trial_id) ; 
				push @phases_drug, $$trial{phase} if ( $d eq Therapy::get_preferred_drug_name( $$trial{'agents'} ) )  ; # && ( $$_{trial_id} ne $trial_id ) 
			}
		}
		my ($max_phases_drug) = sort { $b gt $a } @phases_drug;
		if ( defined $max_phases_drug ) {
			$drug_maturity_tier_cache{$d} = $max_phases_drug ;
			push @phases, $max_phases_drug if defined $max_phases_drug ; # @phases_drug;
		}
	}
	
	my $tier = get_tier_by_phases_of_trials( @phases );
	
	if ($tier =~ /[4]/) {
		my $pubmed_count = 0;
		$pubmed_count += ( $drug_pubmed_count{ Therapy::get_preferred_drug_name( $_ ) } // 0 ) for (@inf_drugs) ;
		$tier .= 'B' if ! $pubmed_count == 0 and $tier !~ /U/;
	}
	return $tier ; 
}

####################################################################################################
my %drug_class_maturity_tier_cache;

sub get_drug_class_maturity_tier {
	&ON_DEMAND_INIT ;
	my $trial_id ;
	$trial_id = shift if ( $_[0] =~ /NCT|ACTRN/ );
	my @inf_drug_classes = @_;
	
	my @phases; 
	
	for my $dc (@inf_drug_classes) {
		next if ! defined $dc ;
		$dc = Therapy::get_normalised_treatment_class_name( $dc );
		if ( exists $drug_class_maturity_tier_cache{$dc} ) {
			push @phases, $drug_class_maturity_tier_cache{$dc};
			next;
		}
		my @phases_drug_classes ;
		if ( grep { Therapy::is_drug_PBS_reimbursed( $_ ) } Therapy::get_all_known_drugs_for_class( $dc ) ) {
			push @phases_drug_classes, 'Phase 4';
			next;
		} elsif ( exists $drug_trial_phases_drug_classes{$dc} ) {
			for my $trial (@{ $drug_trial_phases_drug_classes{$dc} } ) { # @drug_trial_phases_drug_classes
				next if ! defined $$trial{'agents'};
				next if defined $trial_id and ($$trial{'trial_id'} eq $trial_id) ; 
				push @phases_drug_classes, $$trial{phase} if ( $dc eq Therapy::get_normalised_treatment_class_name( $$trial{'agents'} ) )  ; # && ( $$_{trial_id} ne $trial_id ) 
			}
		}
		my ($max_phases_drug_class) = sort { $b gt $a } @phases_drug_classes;
		if ( defined $max_phases_drug_class ) {
			$drug_class_maturity_tier_cache{$dc} = $max_phases_drug_class ;
			push @phases, $max_phases_drug_class if defined $max_phases_drug_class ;
		}
	}
	
	my $tier = get_tier_by_phases_of_trials( @phases );
	
	if ($tier =~ /[4]/) {
		my $pubmed_count = 0;
		$pubmed_count += ( $drug_pubmed_count{ Therapy::get_normalised_treatment_class_name( $_ ) } // 0 )  for (@inf_drug_classes) ;
		$tier .= 'B' if ! $pubmed_count == 0;
	}
	return $tier ; 
}

####################################################################################################
sub get_drug_combination_maturity_tier {
	&ON_DEMAND_INIT ;
	my $trial_id ;
	$trial_id = shift if ( $_[0] =~ /NCT|ACTRN/ );
	my @inf_treatments = @_;
	my @single_agent_tiers;
	my @matched_drug_trial_phases_combo;

	for my $treatment_combo (@inf_treatments) {
		$treatment_combo = Therapy::get_normalised_treatment_name( $treatment_combo ) ;
		next if ! exists $drug_trial_phases_combo{$treatment_combo};
		if ( $treatment_combo !~ /(?: +\+ +|\|)/ ) {;# not a combination
			push @single_agent_tiers, get_drug_maturity_tier( $trial_id, $treatment_combo );
			next;
		}
		
		my @matched_trials = @{ $drug_trial_phases_combo{$treatment_combo} };
		push @matched_trials, @{ $drug_trial_phases_drugs{$treatment_combo} } if exists $drug_trial_phases_drugs{$treatment_combo};
		
		for my $database_combo_row ( @matched_trials ) { # @drug_trial_phases_combo) {
			next if defined $trial_id and ( $$database_combo_row{'trial_id'} eq $trial_id ) ; 
			next if ! defined $$database_combo_row{agents};
			if ( Therapy::get_normalised_treatment_name( $$database_combo_row{agents} ) eq $treatment_combo ) {
				push @matched_drug_trial_phases_combo, $database_combo_row;
				last;
			}
		}
	}
	
# 	print map { join("\t", 'AAA', $$_{'trial_id'}, $$_{'phase'}, $$_{'agents'})."\n" } @matched_drug_trial_phases_combo ;
# 	print map { join("\t", 'BBB', $_, @inf_treatments)."\n" } @single_agent_tiers ;
	
	my ($tier) = sort ( @single_agent_tiers, ( get_tier_by_phases_of_trials( map { $$_{phase} } @matched_drug_trial_phases_combo) ) );
	
	return $tier;
}

####################################################################################################
sub get_drug_class_combination_maturity_tier {
	&ON_DEMAND_INIT ;
	my $trial_id ;
	$trial_id = shift if ( $_[0] =~ /NCT|ACTRN/ );
	my @inf_treatment_class = @_;

	my @single_agent_tiers;
	my @matched_drug_trial_phases_combo_class;
	for my $treatment_combo_class (@inf_treatment_class) {
		$treatment_combo_class = Therapy::get_normalised_treatment_class_name( $treatment_combo_class ) ;
		next if ! exists $drug_trial_phases_combo_classes{$treatment_combo_class} ;
		if ( $treatment_combo_class !~ /(?: +\+ +|\|)/ ) {;# not a combination
			push @single_agent_tiers, get_drug_class_maturity_tier( $trial_id, $treatment_combo_class );
			next;
		}		
		
		my @matched_trials = @{ $drug_trial_phases_combo_classes{$treatment_combo_class} };
		push @matched_trials, @{ $drug_trial_phases_drug_classes{$treatment_combo_class} } if exists $drug_trial_phases_drug_classes{$treatment_combo_class};
		
		for my $database_combo_row ( @matched_trials ) { # @drug_trial_phases_combo_classes 
			next if defined $trial_id and ($$database_combo_row{'trial_id'} eq $trial_id) ; 
			next if ! defined $$database_combo_row{agents};
			if ( Therapy::get_normalised_treatment_class_name( $$database_combo_row{agents} ) eq $treatment_combo_class ) {
				push @matched_drug_trial_phases_combo_class, $database_combo_row;
				last;
			}
		}
	}
	
# 	print map { join("\t", 'CCC', $$_{'trial_id'}, $$_{'phase'}, $$_{'agents'})."\n" } @matched_drug_trial_phases_combo_class ;
# 	print map { join("\t", 'DDD', $_)."\n" } @single_agent_tiers ;

	my ($tier) = sort ( @single_agent_tiers, (get_tier_by_phases_of_trials( map { $$_{phase} } @matched_drug_trial_phases_combo_class ) ) );
	
	return $tier;	
}

####################################################################################################
sub get_drug_maturity_tier_by_trial_id {
	&ON_DEMAND_INIT ;
	my $trial_id = shift;
	my @inf_drugs = @_;
	my @matched = exists $drug_trial_phases_drugs{$trial_id} ? @{ $drug_trial_phases_drugs{$trial_id} } : (); # grep { $$_{trial_id} eq $trial_id } 
	my @tiers ;
	for my $m (@matched) {
		my $flag = 0;
		next if ! length $$m{agents};
# 		print "$$m{trial_id}\t$$m{agents}\n";
		for my $drug (@inf_drugs) {
			if ( Therapy::treatment_contains_drug( $$m{agents}, $drug ) ) {
				$flag = 1 ;
				last;
			}
		}
		next if scalar @inf_drugs and ! $flag;

# 		print join("\t", $$m{trial_id}, $$m{phase}, $$m{agents})."\n";
		push @tiers, get_drug_maturity_tier( $trial_id, $$m{agents} );
	}
	@tiers = sort @tiers;
	return $tiers[0] // '4U';
}


####################################################################################################
sub get_drug_class_maturity_tier_by_trial_id {
	&ON_DEMAND_INIT ;
	my $trial_id = shift;
	my @inf_drugs = @_;
	my @matched = exists $drug_trial_phases_drug_classes{$trial_id} ? @{ $drug_trial_phases_drug_classes{$trial_id} } : (); # grep { $$_{trial_id} eq $trial_id } @drug_trial_phases_drug_classes;
	my @tiers ;
	for my $m (@matched) {
		my $flag = 0;
		next if ! length $$m{agents};
# 		print "$$m{trial_id}\t$$m{agents}\n";
		for my $drug (@inf_drugs) {
			if ( Therapy::treatment_contains_drug( $$m{agents}, $drug ) ) {
				$flag = 1 ;
				last;
			}
		}
		next if scalar @inf_drugs and ! $flag;

# 		print join("\t", $$m{trial_id}, $$m{phase}, $$m{agents})."\n";
		push @tiers, get_drug_class_maturity_tier( $trial_id, $$m{agents} );
	}
	@tiers = sort @tiers;
	return $tiers[0] // '4U';
}

####################################################################################################
sub get_drug_combination_maturity_tier_by_trial_id {
	&ON_DEMAND_INIT ;
	my $trial_id = shift;
	my @inf_drugs = @_;
	my @matched = exists $drug_trial_phases_combo{$trial_id} ? @{ $drug_trial_phases_combo{$trial_id} } : (); # grep { $$_{trial_id} eq $trial_id } @drug_trial_phases_combo;
	
	my @tiers ;
	for my $m (@matched) {
		my $flag = 0;
		next if ! length $$m{agents};
# 		print "$$m{trial_id}\t$$m{agents}\n";
		for my $drug (@inf_drugs) {
			if ( Therapy::treatment_contains_drug( $$m{agents}, $drug ) ) {
				$flag = 1 ;
				last;
			}
		}
		next if scalar @inf_drugs and ! $flag;

# 		print join("\t", $$m{trial_id}, $$m{phase}, $$m{agents})."\n";
		push @tiers, get_drug_combination_maturity_tier( $trial_id, $$m{agents} );
	}
	@tiers = sort @tiers;
	return $tiers[0] // '4U';
}

####################################################################################################
sub get_drug_class_combination_maturity_tier_by_trial_id {
	&ON_DEMAND_INIT ;
	my $trial_id = shift;
	my @inf_drug_class = @_;
	my @matched = exists $drug_trial_phases_combo_classes{$trial_id} ? @{ $drug_trial_phases_combo_classes{$trial_id} } : () ;  # grep { $$_{trial_id} eq $trial_id } @drug_trial_phases_combo_classes;
	my @tiers ;
	for my $m (@matched) {
		my $flag = 0;
		for my $drug_class (@inf_drug_class) {
			if ( Therapy::treatment_contains_drug( $$m{agents}, $drug_class ) ) {
				$flag = 1 ;
				last;
			}
		}
		next if scalar @inf_drug_class and ! $flag;
# 		print join("\t", $$m{trial_id}, $$m{phase}, $$m{agents})."\n";
		push @tiers, get_drug_class_combination_maturity_tier( $trial_id, $$m{agents} );
	}
	@tiers = sort @tiers;
	return $tiers[0] // '4U';
}

######################################################################################################
# assessing maturity of trials 

sub load_drug_trial_phase_database {
	my @dbfiles = @_;
	for my $dbfile (@dbfiles) {
		my @drug_trial_phases = @{ TSV->new($dbfile, [ qw(trial_id phase type agents) ] ) ->{'data'} } ;
		for my $row ( @drug_trial_phases ) {
			next if ! defined $$row{'agents'} or ! length $$row{'agents'} ;
			if ( $$row{'type'} eq 'drug' ) {
				push @{ $drug_trial_phases_drugs{ Therapy::get_preferred_drug_name( $$row{'agents'} ) } }, $row;
				push @{ $drug_trial_phases_drugs{ $$row{'trial_id'} } }, $row;
			} elsif ( $$row{'type'} eq 'drug_class' ) {
				push @{ $drug_trial_phases_drug_classes{ Therapy::get_normalised_treatment_class_name( $$row{'agents'} ) } }, $row;
				push @{ $drug_trial_phases_drug_classes{ $$row{'trial_id'} } }, $row;
			} elsif ( $$row{'type'} eq 'combo' ) {
				push @{ $drug_trial_phases_combo{ Therapy::get_normalised_treatment_name( $$row{'agents'} ) } }, $row;
				push @{ $drug_trial_phases_combo{ $$row{'trial_id'} } }, $row;
			} elsif ( $$row{'type'} eq 'combo_class' ) {
				push @{ $drug_trial_phases_combo_classes{ Therapy::get_normalised_treatment_class_name( $$row{'agents'} ) } }, $row;
				push @{ $drug_trial_phases_combo_classes{ $$row{'trial_id'} } }, $row;
			}
		}
	}
}


###########################################################################################################
sub gen_rules_clinical_trials { # Generating clinical trial rules
	&ON_DEMAND_INIT ;
	my $srcf = shift;
	my $srcf_eligibility = shift; # eligibility file
	my $outf = shift;             # pre-computed cache file
	my $trial_db = TSV->new( $srcf );
	my %trial_db_by_trial_id = $trial_db->index_by('trial_id');
	
	if ( defined($outf) and ( -f $outf ) and ( mtime($srcf) < mtime($outf) ) and 
		( ( defined($srcf_eligibility) and ( mtime($srcf_eligibility) < mtime($outf) ) ) or ! defined($srcf_eligibility) )
		) { 
		return file($outf);
	}
	
	my $Eligibility ;
	my %Eligibility_by_trial_id ;

	if ( defined $srcf_eligibility ) {  # has eligibility file
		$Eligibility = TSV->new( $srcf_eligibility );
		for my $row ( @{ $Eligibility->{'data'} } ) {
			my $trial_id = $$row{'trial_id'};
			my $ec = $$row{'eligibility_criteria'};
			push @{ $Eligibility_by_trial_id{$trial_id} }, $ec;
		}
# 		%Eligibility_by_trial_id = $Eligibility->index_by('trial_id');
	} 

	my @uniq_trials = sort keys %trial_db_by_trial_id;
	my @trials_rules;

	for my $trial_id (@uniq_trials) {
		my $matched_drug_names   = $trial_db_by_trial_id{$trial_id}{'drug_list'}  ;
		my $matched_drug_classes = $trial_db_by_trial_id{$trial_id}{'drug_classes'} ;
		my $trial_acronym        = $trial_db_by_trial_id{$trial_id}{'trialacronym'} // '';
	
		push @trials_rules, mkrule( ['(initial-fact)'],  ["warning:trial $trial_id has no matching drugs."] ) if ! length $matched_drug_names ; 
		
# 		my $healthconditioncode = $trial_db_by_trial_id{$trial_id}{'healthconditioncode'} // '';
# 		my @healthconditioncodes = split /\s*;\s*/, $healthconditioncode ;
		my @healthconditioncodes ;

		my $healthconditions = $trial_db_by_trial_id{$trial_id}{'healthcondition'} // $trial_db_by_trial_id{$trial_id}{'catype'} ;
		my @healthconditions ;
		if ( defined $healthconditions ) {
			my @healthconditions = split /\s*;\s*/, $healthconditions ;
			my @matched_codes = map { DOID::DO_match_catype_whole_word($_) } @healthconditions ;
			@matched_codes = uniq(@matched_codes);
			
# 			print join("\t", $_, DOID::DO_match_catype_whole_word($_) )."\n" for @healthconditions ;
# 			my $s1 = join("\t", sort @matched_codes );
# 			my $s2 = join("\t", sort @healthconditioncodes);
# 			my $clb = ( $s1 eq $s2 ) ? "\e[1;32m" : "\e[1;31m" ;
# 			my $cle = "\e[0m";
# 			
# 			print "HC:  $clb$s1$cle\n";
# 			print "HCD: $clb$s2$cle\n";
			@healthconditioncodes = @matched_codes;
		}
		
		
		my @eligibility_criteria_strs = exists $Eligibility_by_trial_id{$trial_id} ? @{ $Eligibility_by_trial_id{$trial_id} } : ();
		s!(catype:)(.+?)(\s+OR|\s*\)|\s*;|\s*=>|\s*$)!$1.((DO_match_catype($2))[0]//$2).$3!ge for @eligibility_criteria_strs ;
		push @eligibility_criteria_strs, undef if ! scalar @eligibility_criteria_strs;
		
		my @trial_info ;
		push @trial_info, "acronym:".            $trial_acronym ;
		push @trial_info, "drug_list:".          $matched_drug_names ;
		push @trial_info, "drug_classes:".       $matched_drug_classes ;
		push @trial_info, "phase:".              ( $trial_db_by_trial_id{$trial_id}{'phase'}               // '' )  ;
		push @trial_info, "recruitmentstatus:".  ( $trial_db_by_trial_id{$trial_id}{'recruitmentstatus'}   // '' )  ;
		push @trial_info, "approvaldate:".       ( $trial_db_by_trial_id{$trial_id}{'approvaldate'}        // '' )  ;
		push @trial_info, "actualstartdate:".    ( $trial_db_by_trial_id{$trial_id}{'actualstartdate'}     // '' )  ;
		push @trial_info, "studytitle:".         ( $trial_db_by_trial_id{$trial_id}{'studytitle'}          // '' )  ;
		push @trial_info, "healthcondition:".    ( $trial_db_by_trial_id{$trial_id}{'healthcondition'}     // '' )  ;
# 		push @trial_info, "healthconditioncode:".( $trial_db_by_trial_id{$trial_id}{'healthconditioncode'} // '' )  ;
		push @trial_info, "postcode:".           ( $trial_db_by_trial_id{$trial_id}{'postcode'}            // '' )  ;
		push @trial_info, "ext_weblink:".        ( $trial_db_by_trial_id{$trial_id}{'ext_weblink'}         // '' )  ;
		
		my $rhs_str = "preferential_trial_id:$trial_id";
		for my $d ( map { Therapy::get_preferred_drug_name($_) } split /\s*;\s*/, $matched_drug_names )  {
			my @tags = @trial_info ;
			push @tags, "INFERRED:treatment_drug:$d";
			
			if ( scalar(@healthconditioncodes) ) {
				for my $hcc (grep { length } @healthconditioncodes) {
					for my $eligibility_criteria_str (@eligibility_criteria_strs) {
						push @trials_rules, mkrule( 
							[ "sensitive_to:$d", "catype:$hcc",  $eligibility_criteria_str], 
							[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @tags) ]  # 
						) if defined $d;
					}
				}
			} else {
				for my $eligibility_criteria_str (@eligibility_criteria_strs) {
					push @trials_rules, mkrule( 
						[ "sensitive_to:$d", $eligibility_criteria_str ],
						[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @tags) ] # acronym:$trial_acronym; 
					) if defined $d;
				}
			}
		}

		for my $dc ( split /\s*;\s*/, $matched_drug_classes )  {
			for my $dcpart ( split /\s*;\s*/, $dc ) {
				my @tags = @trial_info ;
				push @tags, "INFERRED:treatment_drug_class:$dcpart";
				if ( scalar(@healthconditioncodes) ) {
					for my $hcc (grep { length } @healthconditioncodes) {
						for my $eligibility_criteria_str (@eligibility_criteria_strs) {
							push @trials_rules, mkrule( 
								[ "sensitive_to_drug_class:$dcpart", "catype:$hcc", $eligibility_criteria_str ],
								[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @tags) ] # acronym:$trial_acronym; 
							);
						}
					}
				} else {
					for my $eligibility_criteria_str (@eligibility_criteria_strs) {
						push @trials_rules, mkrule( 
							[ "sensitive_to_drug_class:$dcpart", $eligibility_criteria_str ],
							[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @tags) ] # acronym:$trial_acronym; 
						) ;
					}
				}
			}
		}
	}
	
	@trials_rules = sort( uniq(@trials_rules) );
	
	write_array_uniq($outf, @trials_rules) if defined $outf;

	return @trials_rules;
}

sub load_drug_pubmed_count { 
	my $srcfile = shift;
	for ( file( $srcfile ) ) {
		my ($d, $c, @rest) = split /\t/, $_; 
		$drug_pubmed_count{ Therapy::get_preferred_drug_name($d) } = $c ; 
	}
}



###########################################################################################################
sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;

	for my $srcfile ( POTTRConfig::get_paths('data', 'therapy-maturity-database-file') ) { 
		load_drug_trial_phase_database $srcfile ;
	}
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'therapy-MEDLINE-count-file') ) { 
		load_drug_pubmed_count $srcfile ;
	}
}


######################################################################################################
# SYNOPSIS: 
# 
# print get_drug_maturity_tier('Dabrafenib')."\n";
# 
# print get_drug_class_maturity_tier('BRD_inhibitor')."\n";
# 
# print get_drug_class_maturity_tier('anti-PD-1_monoclonal_antibody', 'PARP_inhibitor')."\n";
# 
# print get_drug_combination_maturity_tier('anti-PD-1_monoclonal_antibody + PARP_inhibitor')."\n";
# 
# print get_drug_combination_maturity_tier('Palbociclib + Fulvestrant')."\n";
# 
# print Therapy::get_normalised_treatment_name('Trametinib + Dabrafenib')."\n";
# 
# print get_drug_class_combination_maturity_tier('anti-PD-1_monoclonal_antibody + PARP_inhibitor')."\n";
# 
# print get_drug_combination_maturity_tier_by_trial_id('NCT02861573', 'Olaparib')."\n";
# 
# print get_drug_class_combination_maturity_tier_by_trial_id('NCT02861573', 'PARP_inhibitor')."\n";
#
# print get_tier_by_phases_of_trials('Phase 2', 'Phase 3', 'Phase 2');
#
# print get_drug_combination_maturity_tier('MK-0482 + Pembrolizumab')."\n";
#
# print get_drug_combination_maturity_tier_by_trial_id('NCT03918278')."\n";
#

1;
