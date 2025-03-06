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


package ClinicalTrials;

use strict;
use warnings;

use POSIX;

use lib '.';
use CancerTypes;
use Evidence;
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
sub uniq { my %order; my %a; my $cnt = 0; for(@_) { $a{$_} = 1; $order{$_} = $cnt if ! exists $order{$_} ; ++$cnt }; return sort {$order{$a} <=> $order{$b}} keys %a; }
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
		m|^(ACTRN\d+[Pp]?)$| and do { return "https://www.anzctr.org.au/TrialSearch.aspx#&&searchTxt=$1" } ;
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


sub clean_catype {
	my $x = $_[0];
	$x =~ s/(?:\b|^)(?:and |or )?((?:any|recurrent|other|relapsed|broad|\s*harbou?ring.*|unspecified|inoperable|locally or metastatic|refractory|adult|childhood|by site|patients? with|cohort.*?:|unresectable|(phase|stage)s? (?:(?:[1234]|I+V?)[ABC]?)|metastatic|MSI-H |locally|advanced|(?:histology|histopathologically) confirmed)(?: and | or |\/| )?)+(?:\b|$)//ig;
	$x =~ s/\s*\([A-Z0-9\-\.]*\)//g;
	$x =~ s/\.$//g;
	return $x;
}



sub gen_rules_catype_match($\@\@) {
	my %retval; 
	my $trial_id = $_[0];
	my @eligibility_criteria_strs = @{ $_[1] };
	my @catypecodes = @{ $_[2] };
	my @catype_general = map { CancerTypes::match_catype_whole_word($_) } ("Neoplasms", "Cancer", "Adenocarcinoma", "Solid tumour", "Liquid cancer", "Haematological malignancy" ) ;
	my $catype_general = join( "|", @catype_general );
	
	my @non_general_hcc = map { /:/ ? $_ : CancerTypes::match_catype_whole_word($_) } grep { length and ( ! /^(?:$catype_general)$/ ) } @catypecodes; 
	
# 	print STDERR "$trial_id: N(eligibility_criteria_strs)=".join(" ", scalar(@eligibility_criteria_strs) )."\n";
	for my $eligibility_criteria_str (@eligibility_criteria_strs) {
	
# 		next if ! length $eligibility_criteria_str;
		
		my $f_eligibility_criteria_str_has_catype = ( (length $eligibility_criteria_str) && ( $eligibility_criteria_str =~ /(?:^|;\s*|\bOR\s*)catype:/ ) );
				
		$retval{$eligibility_criteria_str} = [] if ! exists $retval{$eligibility_criteria_str} ;
		
		my @eligibility_criteria_str_catypes = map { /^\(.*\bor\b.*\)$/i ? ( split /\s+OR\s+/i, ( s/\bcatype:|^\(|\)$//gr ) ) : (s/^catype://r) } grep { /^\(?catype:/ } split /\s*;\s*/, ($eligibility_criteria_str // '');
# 		print STDERR "\t>$_\n" for @eligibility_criteria_str_catypes ;
				
		my $cnt = 0;
		for my $hcc (grep { $f_eligibility_criteria_str_has_catype } @non_general_hcc) {
# 			print STDERR ">> $trial_id || $f_eligibility_criteria_str_has_catype \t|| [$hcc] || ".($eligibility_criteria_str//'')." \n";
					
			for my $ecc ( @eligibility_criteria_str_catypes ) {
				my $hcc_isa_ecc = CancerTypes::is_a( $hcc, $ecc );
				my $ecc_isa_hcc = CancerTypes::is_a( $ecc, $hcc );
				my $flag = ( $hcc_isa_ecc * $ecc_isa_hcc ) ? "=" : ( $hcc_isa_ecc ? "(" : ( $ecc_isa_hcc ? ')' : "X" ) );
				my $col  = ( $hcc_isa_ecc * $ecc_isa_hcc ) ?  33 : ( $hcc_isa_ecc ? 36 : ( $ecc_isa_hcc ? '38;5;72' : 30) );
# 				print STDERR "\t\e[1;${col}m$hcc $flag $ecc\e[0m\n";
				++$cnt;
				next if $flag eq 'X'; # $flag eq ')' 
				if ( $flag eq '=' or $flag eq ')' ) {
					push @{ $retval{$eligibility_criteria_str} }, '';
				} elsif ( $flag eq '(' ) {
					push @{ $retval{$eligibility_criteria_str} }, join("; ", "catype:$ecc", "*catype:$hcc") if $ecc ne 'Cancer';
				}
			}
		}
		
		if ( ! $cnt ) {
# 			print STDERR "\t> !\$cnt\tnon_general_hcc=".join(", ", @non_general_hcc)."\n" ;
		
			for my $hcc (@non_general_hcc) {
				push @{ $retval{$eligibility_criteria_str} }, "catype:$hcc" ;
# 				print STDERR "\t> added $hcc\n" ;
				my @hcc_ancestors = ( $hcc =~ /:/ ? ($hcc) : ( grep { defined and ! /Solid|^Cancer$|Haematologic/ } CancerTypes::get_ancestors( $hcc ) ) );
				my @hcc_ancestor_is_a_reference_catype = grep { defined and exists $Evidence::reference_catypes{$_} } @hcc_ancestors;
				if ( scalar(@hcc_ancestor_is_a_reference_catype) ) {
					A: for my $a (@hcc_ancestors) {
# 						print STDERR "$trial_id\t$hcc\t$a\t".scalar(@hcc_ancestor_is_a_reference_catype)." ".join(" ", @hcc_ancestor_is_a_reference_catype)."\n";
						next A if $a eq $hcc;
						for my $b (@hcc_ancestor_is_a_reference_catype) {
							my $a_isa_b = CancerTypes::is_a( $a, $b );
							next A if $a_isa_b ;
						}
						push @{ $retval{''} }, "catype:$a; *catype:$hcc" ;
					}
				}
			}
		}
		
		push @{ $retval{$eligibility_criteria_str} }, $eligibility_criteria_str if ( $eligibility_criteria_str =~ /catype:/ ) and ( ! grep { $_ eq $eligibility_criteria_str } @{ $retval{$eligibility_criteria_str} } );
	}
	return %retval;
}

sub gen_rules_clinical_trials($\@$) { # Generating clinical trial rules
	&ON_DEMAND_INIT ;
	my $srcf = $_[0];
	my @srcf_eligibility = @{ $_[1] } ; # eligibility files
	my $outf = $_[2];                   # pre-computed cache file for the rules
	my $trial_db = TSV->new( $srcf );
	my %trial_db_by_trial_id = $trial_db->index_by('trial_id');
	
	if ( defined($outf) and ( -f $outf ) and ( mtime($srcf) < mtime($outf) ) ) { 
		my $cache_is_recent = 1;
		for my $srcf_eligibility ( @srcf_eligibility ) {
			next if ! defined($srcf_eligibility) ;
			next if ! -f $srcf_eligibility;
			next if mtime($srcf_eligibility) < mtime($outf) ;
			$cache_is_recent = 0;
			last;
		}
		return file($outf) if $cache_is_recent ;
	}
	
	my %Eligibility_by_trial_id ;
	
	# The final eligibility of a trial is populated in conjunction of rule sets defined in each file (whereas 
	# within each rule set, the rules are populated disjunction). Specifically: 
	#   eligibility_list_1 x eligibility_list_2 x ... x eligibility_list_n =
	#   { e11 & e21, e11 & e22, ... e12 & e21, e212 &e22 ...  }
	# where eligibility_list_1 = { e11, e12, ... e1k }, eligibility_list_2 = { e21, e22, ... e2m } 
	
	@srcf_eligibility = uniq(@srcf_eligibility);
	
# 	die join("\n", @srcf_eligibility)."\n";
	
	for my $srcf_eligibility ( @srcf_eligibility ) {
		my $Eligibility_list_in_this_file = ( ($srcf_eligibility =~ /.tsv/) ? 
			TSV->new( $srcf_eligibility ) : # if defined as a TSV, then load as is, 
			TSV->new( $srcf_eligibility, ['trial_id', 'eligibility_criteria'] ) # otherwise assume default order of fields of trial_id, ec
		);
			
		my %Eligibility_by_trial_id_this_file ;
		for my $row ( @{ $Eligibility_list_in_this_file->{'data'} } ) {
			my $trial_id = $$row{'trial_id'};
			my $ec = $$row{'eligibility_criteria'};
			next if ! defined $trial_id;
			next if ! defined $ec;
			# FIXME: A temporarily fix to whitelist the recruitment criteria
			next if $ec =~ /^recruiting\s*(?:$|\b|\[)/i;
			push @{ $Eligibility_by_trial_id_this_file{$trial_id} }, $ec; # rules are populated in disjunction
		}
		
		for my $trial_id ( keys %Eligibility_by_trial_id_this_file ) {
			if ( exists $Eligibility_by_trial_id{$trial_id} ) { # if old rule set already exists, concatenate the new rules behind each of the old rules
				for my $ec_old ( @{ $Eligibility_by_trial_id{$trial_id} } ) {
					for my $ec_new ( @{ $Eligibility_by_trial_id_this_file{$trial_id} } ) {
						$ec_old = join('; ', (  ( $ec_old, $ec_new) ) ); # grep { defined }
					}
				}
			} else { # populate all rules in disjunction.
				push @{ $Eligibility_by_trial_id{$trial_id} }, @{ $Eligibility_by_trial_id_this_file{$trial_id} }; 
			}
		}
	}

	my @uniq_trials = sort keys %trial_db_by_trial_id;
	my @trials_rules;

	for my $trial_id (@uniq_trials) {
		my $trial_acronym = $trial_db_by_trial_id{$trial_id}{'trialacronym'} // '';
		
		my $matched_drug_names    = $trial_db_by_trial_id{$trial_id}{'drug_list'}  ;
		my $matched_drug_classes  = $trial_db_by_trial_id{$trial_id}{'drug_classes'} ;
		my $matched_combos        = $trial_db_by_trial_id{$trial_id}{'combo_list'}  || $trial_db_by_trial_id{$trial_id}{'drug_list'}  ;
		my $matched_combo_classes = $trial_db_by_trial_id{$trial_id}{'combo_classes'} || $trial_db_by_trial_id{$trial_id}{'drug_classes'} ;
		
		my @matched_drug_names = map { Therapy::get_normalised_treatment_name($_) } ( grep { length } split /\s*;\s*/, $matched_drug_names ) ;
		my @matched_drug_classes = split /\s*;\s*/, $matched_drug_classes ;
		my @matched_combos = map { Therapy::get_normalised_treatment_name($_) } ( grep { length } split /\s*;\s*/, $matched_combos ) ;
		my @matched_combo_classes = split /\s*;\s*/, $matched_combo_classes ;
		
# 		push @trials_rules, mkrule( ['(initial-fact)'],  ["warning:trial $trial_id has no matching drugs."] ) if ! scalar @matched_drug_names; 

		my @catypecodes ;

		my $healthconditions = $trial_db_by_trial_id{$trial_id}{'healthcondition'} // $trial_db_by_trial_id{$trial_id}{'catype'} ;
		my @healthconditions ;
		if ( defined $healthconditions ) {
			my @healthconditions = split /\s*;\s*/, $healthconditions ;

			for (@healthconditions) {
				s/(?:metastatic|advanced)\s+//i ;
				s/, adult$//gi;
				s/tumor/tumour/gi;
			}
			
			my @matched_catypes = map { 
				my $mm ;
				my $retval =
					/^(\S+)[ \-]+(?:Gene )?(Mutations|mutant)?$/i  ?  "$1:oncogenic_mutation" :
					/^(\S+)[ \-]+(?:Gene )?Alterations?$/i ?  "$1:alteration" :
					/^(\S+)[ \-]+(?:Gene )?Fusions?$/i ?  "$1:fusion" :
					/^(\S+)[ \-]+(?:Wild.?type)/i ?  "NOT $1:oncogenic_mutation" :
					/^(\S+)[ \-]+(?:Gene )?(?:deletion|deleted)?$/i ?  "$1:deletion" :
					($mm = CancerTypes::match_catype_whole_word( clean_catype($_), 'STRICT' )) ;
					
				my @xx;
				if ( ! $mm ) {
					my %xx = CancerTypes::match_catype($_);
					@xx = keys %xx;
# 					print STDERR "\e[1;35mCancerTypes.pm: \e[0m$_ => \e[1;35m".join("; ", @xx)." \e[0m\n";
				}
				$retval // @xx;
				} @healthconditions;
				
			@catypecodes = uniq(@matched_catypes);
			
# 			print join("\t", "T", $trial_id, $_, CancerTypes::match_catype_whole_word($_) )."\n" for @healthconditions ;
# 			my $s1 = join("\t", sort @matched_codes );
# 			my $s2 = join("\t", sort @catypecodes);
# 			my $clb = ( $s1 eq $s2 ) ? "\e[1;32m" : "\e[1;31m" ;
# 			my $cle = "\e[0m";
# 			
# 			print "\tHC:  $clb$s1$cle\n";
# 			print "\tHCD: $clb$s2$cle\n";
		}
# 		print join("\t", "T", $trial_id, @catypecodes )."\n" ;
		
		my @eligibility_criteria_strs = exists $Eligibility_by_trial_id{$trial_id} ? @{ $Eligibility_by_trial_id{$trial_id} } : ();
		s!(catype:)(.+?)(\s+OR|\s*\)|\s*;|\s*=>|\s*$)!$1.((CancerTypes::match_catype($2))[0]//$2).$3!ge for @eligibility_criteria_strs ;
		push @eligibility_criteria_strs, '' if ! scalar @eligibility_criteria_strs;  # add an undefined value if empty eligibility string
		
		my @trial_info ;
		push @trial_info, "acronym:".            $trial_acronym ;
		push @trial_info, "drug_list:".          $matched_drug_names ;
		push @trial_info, "drug_classes:".       $matched_drug_classes ;
		push @trial_info, "combo_list:".         $matched_combos ;
		push @trial_info, "combo_classes:".      $matched_combo_classes ;
		push @trial_info, "phase:".              ( $trial_db_by_trial_id{$trial_id}{'phase'}               // '' )  ;
		push @trial_info, "recruitmentstatus:".  ( $trial_db_by_trial_id{$trial_id}{'recruitmentstatus'}   // '' )  ;
		push @trial_info, "approvaldate:".       ( $trial_db_by_trial_id{$trial_id}{'approvaldate'}        // '' )  ;
		push @trial_info, "actualstartdate:".    ( $trial_db_by_trial_id{$trial_id}{'actualstartdate'}     // '' )  ;
		push @trial_info, "studytitle:".         ( $trial_db_by_trial_id{$trial_id}{'studytitle'}          // '' )  ;
		push @trial_info, "healthcondition:".    ( $trial_db_by_trial_id{$trial_id}{'healthcondition'}     // '' )  ;
		push @trial_info, "postcode:".           ( $trial_db_by_trial_id{$trial_id}{'postcode'}            // '' )  ;
		push @trial_info, "ext_weblink:".        ( $trial_db_by_trial_id{$trial_id}{'ext_weblink'}         // '' )  ;

		my $rhs_str = "preferential_trial_id:$trial_id";

		push @trials_rules, mkrule( [ $rhs_str ], [ Facts::mk_fact_str( $rhs_str, @trial_info) ] ) ; # autoannotating trials
		
		my $n_preferential_rules = 0; # number of preferential trials listed
		
		my $catypecodes_str = ( scalar(@catypecodes) ? join(' OR ', ( map { /mutation/i ? $_ : "catype:$_" } @catypecodes )) : undef ) ;
		$catypecodes_str = "($catypecodes_str)" if defined $catypecodes_str and $catypecodes_str =~ / OR /;
		
		 # adding information about recruitment status or available - consistent with CT.gov terminology
		if ( $trial_db_by_trial_id{$trial_id}{'recruitmentstatus'} !~ /^(?:Recruiting|Available)$/i ) {
			if ( scalar(@eligibility_criteria_strs) ) {
 				$_ .= ( length $_ ? "; " : "" )."*recruitment_status:".$trial_db_by_trial_id{$trial_id}{'recruitmentstatus'} for @eligibility_criteria_strs;
				; 
			} else {
				push @eligibility_criteria_strs, "*recruitment_status:".$trial_db_by_trial_id{$trial_id}{'recruitmentstatus'}; 
			}
		}
		
		my @tags = @trial_info ;

		my %extra_catype_from_eligibility_str = gen_rules_catype_match($trial_id, @eligibility_criteria_strs, @catypecodes); 
# 		print STDERR "| $trial_id\n";

		for my $eligibility_criteria_str (@eligibility_criteria_strs) { # for each manually defined eligibility criteria (empty string if none)
		
			my $f_eligibility_criteria_str_has_sensitivity = ( (length $eligibility_criteria_str) && ( $eligibility_criteria_str =~ /(?:^|;\s*|\bOR\s*)sensitive_to(?:_drug_class)?:/ ) );
			my $f_eligibility_criteria_str_has_catype      = ( (length $eligibility_criteria_str) && ( $eligibility_criteria_str =~ /(?:^\(?\s*|;\s*|\bOR\s*)catype:/ ) );

			my @extra_catype_strs = length $eligibility_criteria_str ? @{ $extra_catype_from_eligibility_str{$eligibility_criteria_str} } : ();
			
# 			print "|| $trial_id | f_eligibility_criteria_str_has_sensitivity=$f_eligibility_criteria_str_has_sensitivity\teligibility_criteria_str=$eligibility_criteria_str\n";
# 			print "|| >>> extra_catype_strs=[".join(" ", @extra_catype_strs)."] \n";
			
			for my $d ( @matched_drug_names )  {
				my $trial_match_criteria_str = "trial_match_criteria:drug_sensitivity";
				for my $c ( @matched_combos )  {
# 					print STDERR "|| $d\t$c\n";
					my $drug_sensitivity_str = "sensitive_to:$d; NOT therapy_recommendation:$d:U";
					my $drug_sensitivity_tag = "INFERRED:treatment_drug:$d";
					my $f_preferential_trial_combo_complete_match = undef;
					if ( $c eq $d ) {
						$f_preferential_trial_combo_complete_match = 'preferential_trial_complete_match:sensitive_to';
					} elsif ( $c =~ /(?:^|\b)\Q$d\E(?:$|\b)/ ) {
						$f_preferential_trial_combo_complete_match = 'preferential_trial_partial_match:sensitive_to';
						$drug_sensitivity_str .= "; *sensitive_to:$c";
					} else {
						next;
					}

					my $f_eligibility_criteria_str_has_drug_sensitivity = ( $f_eligibility_criteria_str_has_sensitivity and ( $eligibility_criteria_str =~ /(?:^|;\s*|\bOR\s*)$drug_sensitivity_str(?:\d*$|\b)/ ) );
					my $f_assert_drug_sensitivity_tag = ( (! $f_eligibility_criteria_str_has_sensitivity) or $f_eligibility_criteria_str_has_drug_sensitivity) ;
				
# 					print STDERR "||| $f_eligibility_criteria_str_has_sensitivity\t$drug_sensitivity_str\t|\t$eligibility_criteria_str\n";

					for my $extra_catype_strs ('', @extra_catype_strs) {
						next if $extra_catype_strs eq ($catypecodes_str // '');
						my $f_preferential_trial_catype_complete_match = ( join('; ', $extra_catype_strs, ($catypecodes_str // ''), $eligibility_criteria_str // '') !~ /\*catype:/ ? 'preferential_trial_complete_match:catype' : undef);
# 						my $trial_match_criteria_str_spec = $trial_match_criteria_str;
# 						$trial_match_criteria_str_spec = join("; ", $trial_match_criteria_str_spec, "trial_match_criteria:cancer_type") if length $extra_catype_strs;
						my @trial_match_criteria_str_spec = ($trial_match_criteria_str);
						push @trial_match_criteria_str_spec, "trial_match_criteria:cancer_type" if ( length $extra_catype_strs || $f_eligibility_criteria_str_has_catype ); 

						push @trials_rules, mkrule( 
							[ ( $f_eligibility_criteria_str_has_sensitivity ? undef : $drug_sensitivity_str) , 
								( $f_eligibility_criteria_str_has_catype ? undef : ( ( $extra_catype_strs =~ /\*$catypecodes_str(?:;|$)/) ? undef : $catypecodes_str )  ) ,
								$extra_catype_strs, 
								$eligibility_criteria_str
							], 
							[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @trial_info, @trial_match_criteria_str_spec, 'drug_specific_rule',
								($f_assert_drug_sensitivity_tag ? $drug_sensitivity_tag: undef),
								$f_preferential_trial_catype_complete_match,
								$f_preferential_trial_combo_complete_match 
								),
							]  # 
						) ;
					}

# 					print "|||>  $trials_rules[$#trials_rules]\n";
					++ $n_preferential_rules;
				}
			}
			
# 			print STDERR "| DC: ".join(" ", @matched_drug_classes)."\n";
			for my $dc ( @matched_drug_classes )  {
				my @tags = @trial_info ;
				for my $cc ( @matched_combo_classes )  {
# 					print STDERR "| DC: CC: $cc\n";
					my $trial_match_criteria_str = "trial_match_criteria:drug_class_sensitivity";
						
					my $drug_class_sensitivity_str = "sensitive_to_drug_class:$dc";
					my $drug_class_sensitivity_tag = "INFERRED:treatment_drug_class:$dc";

					my $f_preferential_trial_combo_class_complete_match = undef;
					if ( $cc eq $dc ) {
						$f_preferential_trial_combo_class_complete_match = 'preferential_trial_complete_match:sensitive_to_drug_class';
					} elsif ( $cc =~ /(?:^|\b)\Q$dc\E(?:$|\b)/ ) {
						$drug_class_sensitivity_str .= "; *sensitive_to_drug_class:$cc";
					} else {
						next;
					}
					
					
					my $f_eligibility_criteria_str_has_drug_class_sensitivity = ( $f_eligibility_criteria_str_has_sensitivity and ( $eligibility_criteria_str =~ /(?:^|;\s*|\bOR\s*)$drug_class_sensitivity_str(?:\d*$|\b)/ ) );
					my $f_assert_drug_class_sensitivity_tag = ( (! $f_eligibility_criteria_str_has_sensitivity) or $f_eligibility_criteria_str_has_drug_class_sensitivity) ;
							
# 					print STDERR "||| DC: f_eligibility_criteria_str_has_drug_class_sensitivity=$f_eligibility_criteria_str_has_drug_class_sensitivity\tdrug_class_sensitivity_str=$drug_class_sensitivity_str\t|\teligibility_criteria_str=$eligibility_criteria_str\n";

					for my $extra_catype_strs ('', @extra_catype_strs) {
						next if $extra_catype_strs eq ($catypecodes_str // '');
						my $f_preferential_trial_catype_complete_match = ( join('; ', $extra_catype_strs, ($catypecodes_str // ''), $eligibility_criteria_str // '') !~ /\*catype:/ ? 'preferential_trial_complete_match:catype' : undef);
						my @trial_match_criteria_str_spec = ($trial_match_criteria_str);
						push @trial_match_criteria_str_spec, "trial_match_criteria:cancer_type" if ( length $extra_catype_strs || $f_eligibility_criteria_str_has_catype ); 
						
						push @trials_rules, mkrule( 
							[ ( $f_eligibility_criteria_str_has_sensitivity ? undef : $drug_class_sensitivity_str), 
								( $f_eligibility_criteria_str_has_catype ? undef : ( ( $extra_catype_strs =~ /\*$catypecodes_str(?:;|$)/) ? undef : $catypecodes_str ) ), 
								$extra_catype_strs,
								$eligibility_criteria_str 
							],
							[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @trial_info, @trial_match_criteria_str_spec, 'drug_class_specific_rule',
								($f_assert_drug_class_sensitivity_tag ? $drug_class_sensitivity_tag : undef),
								$f_preferential_trial_catype_complete_match,
								$f_preferential_trial_combo_class_complete_match
								)
							] # acronym:$trial_acronym; 
						);
					}
					++ $n_preferential_rules;
				}
			}
		}


# 		print "$n_preferential_rules\t$matched_drug_names\n";
# 		print STDERR "$trial_id:\t".join(" ", @catypecodes)."\n";
		if ( scalar(@catypecodes) ) {
			my @tags = @trial_info ;
			push @tags, "trial_match_criteria:cancer_type";
			for my $eligibility_criteria_str (@eligibility_criteria_strs) {
				my @extra_catype_strs = @{ $extra_catype_from_eligibility_str{$eligibility_criteria_str} } ;
				
# 				print STDERR "$trial_id:\teligibility_criteria_str=$eligibility_criteria_str: ".join(" ", scalar(@eligibility_criteria_strs), scalar(@extra_catype_strs) )."\n";
			
# 				print STDERR "\t\e[1;38;5;97mhealthconditions      =$healthconditions\e[0m\n";
# 				print STDERR "\t\e[1;38;5;97m    extra_catype_strs =$_\e[0m\n" for @extra_catype_strs ;
				for my $extra_catype_str (@extra_catype_strs) {
# 					print STDERR "\t\e[1;38;5;144mextra_catype_str=$extra_catype_str\e[0m; \e[1;38;5;64meligibility_criteria_str=$eligibility_criteria_str\e[0m\n";
					push @trials_rules, mkrule( 
						[ $extra_catype_str, $eligibility_criteria_str ], 
						[ Facts::mk_fact_str("preferential_trial_id:$trial_id", @tags, 'cancer_type_specific_rule',
							( join('; ', $extra_catype_str, $eligibility_criteria_str // '') !~ /\*catype:/ ? 'preferential_trial_complete_match:catype' : undef)
							)
						]  # 
					) ;
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
