#!/usr/bin/perl
##############################################################################
#
# pottr.pl: Precision Oncology Therapy and Trial Recommender -
# 
# A decision support system for guiding molecularly targeted cancer therapies
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
#
# Version history: 
#   Version 0.5  - 22 May 2019 
#   Version 0.5  - 30 Sep 2019 
#   Version 0.6  - 18 Nov 2019 
#   Version 0.7  - 26 Nov 2019 
#   Version 0.9  - 20 May 2020
#   Version 0.91 - 27 May 2020
#   Version 0.92 - 13 Jun 2020
#   Version 0.93 - 14 Sep 2020
#   Version 0.95 - 24 Dec 2021
#
##############################################################################

our $POTTR_version = '0.95';
our $POTTR_title = qq|Precision Oncology Treatment and Trial Recommender (Version $POTTR_version)|;

our $POTTR_subtitle = qq|
A clinical decision support system for guiding targeted therapies in cancer

POTTR on github: http://github.org/fpylin/POTTR
|;

our $POTTR_disclaimer = qq|
DISCLAIMER

Specialist knowledge in oncology is constantly evolving. Making treatment recommendations  
for cancer patients requires  careful evaluation  of clinical situations in its entirety.  
Results produced by POTTR is  STRICTLY FOR RESEARCH USE ONLY  and  CANNNOT SUBSITUTE  for 
professional oncology advice. Using POTTR is WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND. 
Inappropriate use of this software may lead to harm.

This is a beta software.  For bug fixes, please contact: Dr. Frank Lin, MBChB, PhD, FRACP
(f dot lin at garvan dot org dot au).
|;

use strict;
use warnings;


use Getopt::Long qw(:config bundling);
use Pod::Usage;

use File::Basename;
use lib dirname (__FILE__);

use POTTRConfig ( dirname (__FILE__) );

use TSV;
use Rules;
use POTTR;
use POSIX;
use Data::Dumper;

#######################################################################

sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub max { return undef if (! scalar(@_) ); my $v = shift; for (@_) { $v = $_ if ($_ > $v) ; } return $v; }
sub column {
	my @lines = @_;
	my @levels;
	my @width;
	
	my $retval;
	for (@lines) {
		chomp;
		my @parts = split /\t/, $_;
		$levels[$_]{ $parts[$_] }++ for (0 .. $#parts);
		}
		
	sub deansi { my $x = shift; $x =~ s/(\e\[[0-9;]+m)//g ; return $x;}
	
	for my $i (0 .. $#levels) {
# 		$width[$i] = max( map { my $a = $_; $a =~ s/\e\[[0-9;]+m//g; length($a) } keys %{ $levels[$i] } );
		$width[$i] = max( map { my $a = $_; length(deansi($a)) } keys %{ $levels[$i] } );
# 		$width[$i] = max( map { my $a = $_; length($a) } keys %{ $levels[$i] } );
		}
	

	sub normANSIcc {
		my $x = shift;
		my $max_width = shift;
		my $append = '';
		my $cc_length = 0;
		
		while ( $x =~ /(\e\[[0-9;]+m)/g ) {
			my $ctrl_code = $1;
			$cc_length += length($ctrl_code);
			}
			
			my $deansied_x = deansi($x); 
			if ( (length($deansied_x) + $cc_length) > $max_width ) {
				$cc_length  = $max_width - length($deansied_x);
				}
		
			$append .= " " x $cc_length if $cc_length > 0;
			return $x.$append;
	}

	sub adj {
		my $width = shift;
		my $str = shift;
		$str = normANSIcc($str, $width);
		$str .= " " x max($width - length(deansi($str)), 0);
		return $str;
		}
	
	for (@lines) {
		my @parts = split /\t/, $_;
		my @output = map { adj($width[$_], $parts[$_]) } (0..$#parts);
		my $line = join("   ", @output);
		$retval .= $line."\n";
		}
	return $retval;
	}


sub print_column {
	my $x = shift;
	our @buffer;
	if ( ! defined $x ) {
		print for column @buffer;
		@buffer=();
	} else {
		push @buffer, $x;
	}
}
	

sub column_wrap_text($@) {
	my $width = shift;
	my $width_minus_5 = $width - 5;
	my @lines = @_;
	
	@lines = map { # Wrapping function
		chomp;
		my @in = split /\t/, $_;
		my @out;
		my $r = 0;
		for my $c ( 0 .. $#in ) {
			my $skip = 0;
			my $x = $in[$c];
			$x =~ s/\s*$//;
	#      		print ">>[$x]\n";
			while ( length($x) >= $width ) {
				my ($a, $b) = ( $x =~ /^(.{$width_minus_5,}?[^[:space:];]+[\s;]*)(.*)/ );
	# 				print "$x -> [$a] [$b]\n";
				if ( $b =~ /^\s*$/ ) { $x = $a; last; }
				$x = $b;
				next if $a =~ /^\s*$/;
				$out[$r+$skip][$c] = $a;
	# 				last if $b =~ /^\s*$/;
				++$skip;
			} 
	# 			print "  [$x]\n";
			$out[$r+$skip][$c] = $x ; # if $x !~ /^\s*$/;
			$r += $skip ; # + 1;
		}
		( map { join("\t", ( map { ($_ // '') } @{$_} ) ) } @out );
		} @lines ;
	return @lines;
}
#######################################################################



#######################################################################
# my $f_from_web = ( (exists $ENV{'HOME'}) ) ? 0 : 1 ;
my $f_from_web = 0 ;

my $pottr_modules = undef;
my $f_webpage = undef;
my $f_title = undef;
my $f_quiet = 1;
my $f_verbose = undef;
my $f_print_rules = undef;
my $f_tsv = undef;

my $f_config = "pottr_config";

my $f_interp_catype = undef;
my $f_interp_variants = undef;
my $f_list_trials = undef;
my $f_list_evidence = undef;
my $f_list_therapies = undef;
my $f_list_therapy_summaries = undef;
my $f_clinical_report = undef;
my $f_help = undef;
my $f_man = undef;
my $f_version = undef;

my $f_list_trials_terse = undef; # print less information about trials ;

GetOptions (
	"help|h|?"                => \$f_help,
	"man|manual|docs|m"       => \$f_man,
	"web-output|w"            => \$f_from_web,
	"webpage"                 => \$f_webpage,
	"title=s"                 => \$f_title,
	"modules=s"               => \$pottr_modules,
	"quiet|q"                 => \$f_quiet,
	"config-file|c=s"         => \$f_config,
	"version"                 => \$f_version,
	"verbose|debug|d|v"       => \$f_verbose,
	"print-rules|r"           => \$f_print_rules,
	"clinical-report|p"       => \$f_clinical_report,
	"export-tsv|t"            => \$f_tsv,
	"interpret-catypes|C"     => \$f_interp_catype,
	"interpret-variants|V"    => \$f_interp_variants,
	"list-evidence|E"         => \$f_list_evidence,
	"list-therapies|R"        => \$f_list_therapies,
	"list-therapy-summary|S"  => \$f_list_therapy_summaries,
	"list-trials|T"           => \$f_list_trials,
	"trials-terse"            => \$f_list_trials_terse,
);

pod2usage(-verbose=>1) if $f_help;
pod2usage(-verbose=>99) if $f_man;
if ($f_version) {
	print "Precision Oncology Treatment and Trial Recommender $POTTR_version\n";
	exit(0);
}


if ( ! ($f_interp_catype || $f_interp_variants || $f_list_evidence || $f_list_trials || $f_list_therapy_summaries|| $f_list_therapies) ) {
	$f_interp_catype = $f_interp_variants = $f_list_evidence = $f_list_therapy_summaries = $f_list_trials = $f_list_therapies = 1;
}

POTTRConfig::load( $f_config ); # loading configuration file

$f_quiet = 0 if $f_verbose ;

my $Pottrparams = {
	'debug_level' => ( ( ($f_quiet) || $f_from_web ) ? 0 : 1  ), 
	'debug_msg'   => ( $f_from_web ? 'html' : 'ansi' )
	} ;
	
my @pottr_modules;
if (! defined $pottr_modules) {
	push @pottr_modules, "CTM" if defined $f_interp_catype ;
	push @pottr_modules, "VFM" if defined $f_interp_variants ;
	push @pottr_modules, "CTM CLI VFM VEG" if defined $f_list_evidence ;
	push @pottr_modules, "CTM CLI VFM VEG DSO" if defined $f_list_therapies ;
	push @pottr_modules, "CTM CLI VFM VEG DSO TRG" if defined $f_list_therapy_summaries ;
	push @pottr_modules, "CTM CLI VFM VEG DSO TRG PTM PTP" if defined $f_list_trials ;
}

$Pottrparams->{'modules'} = ( (defined $pottr_modules) ? $pottr_modules : join(" ", @pottr_modules) );

my $Pottr = POTTR->new( $Pottrparams );

if ($f_print_rules) {
	$Pottr->{'modules'}->print();
	exit 0;
}

############################################################################

my @lines = scalar @ARGV ? @ARGV : <STDIN>;
my @input = @lines ;
s/^\s+|\s+$//g for @lines ;
@lines = grep { length $_ } @lines;
Biomarker::load_biomarker_synonym_file;
for (@lines) {
	chomp;
	next if /:/;
	$_ = "$1:$2" if ( /^(\S+)?\s+(.+)/ and ( exists $Biomarker::biomarker_synonyms{uc($1)} or $1 eq 'catype') );
}
# unshift @lines, "catype:Solid tumour" if ! grep {/catype:/} @lines;
my @facts = grep { length } @lines;
my %facts = map { s/\s*\[.*//r => 1 } @facts;

my @treatment_match_result = $Pottr->reason(@facts);

# exit 0  if (! $f_from_web and ! $export_trial_list_filename );

#######################################################################
my @results_catypes              = grep { /catype_name:/ }                                         @treatment_match_result ;
my @results_biomarker            = grep { /has_biomarker:/ }                                       @treatment_match_result ;
my @results_treatment            = grep { /(?:^\S+:treatment:)/ }                                  @treatment_match_result ;
my @results_sens_to_drugs        = grep { /^(?:recommendation_tier:[^:]+?:[^R])/ }                 @treatment_match_result ;
my @results_resi_to_drugs        = grep { /^(?:recommendation_tier:[^:]+?:R)/ }                    @treatment_match_result ;
my @results_sens_drug_class      = grep { /^(?:recommendation_tier_drug_class:[^:]+?:[^R])/ }      @treatment_match_result ;
my @results_resi_drug_class      = grep { /^(?:recommendation_tier_drug_class:[^:]+?:R)/ }         @treatment_match_result ;
my @results_therapy_recommendations = grep { /^(?:therapy_recommendation:)/ }                          @treatment_match_result ;
my @results_preferential_trials  = grep { /^(?:preferential_trial_id:)/ }                          @treatment_match_result ;

#######################################################################################################################################
sub extract_LOE {
	return $1 if $_[0] =~ /(?:LOE: |^recommendation_tier\S*?:.*?:)(R?[01234S][AB]?)/;
	return  0 ;
}

sub extract_LOM {
	return $1 if $_[0] =~ /LOM:([123])/;
	return 9;
}

sub extract_pref_score {
	if ( $_[0] =~ /(?:^recommendation_tier\S*?:.*?:).*\bpref_score:(\d+(?:\.\d+)?)/ ) {
		return $1 
	}
	return  0 ;
}

sub score_tier { return POTTR::score_tier(extract_LOE($_[0])); }

sub score_LOM { return extract_LOM($_[0]); }

sub hl { my ($col, $x) = @_; return "\e[1;${col}m$x\e[0m"; }

sub thl {
	my $tier = shift;
	my $col  = 37;
	$col = 32 if $tier =~ /^[1]/;
	$col = 36 if $tier =~ /^[2]/;
	$col = 33 if $tier =~ /^[3]/;
	$col = 35 if $tier =~ /^[4]/;
	$col = 31 if $tier =~ /^R/;
	$col = 30 if $tier =~ /U/;
	$col = 30 if $tier =~ /(?:^|,|\b)(?:\(|--)/;
	return hl($col, $tier);
}

sub thl_AMP {
	my $tier = shift;
	my $col  = 37;
	$col = 32 if $tier =~ /^1A/;
	$col = 36 if $tier =~ /^1B/;
	$col = 35 if $tier =~ /^2C/;
	$col = 34 if $tier =~ /^2D/;
	$col = 30 if $tier =~ /^[34]/;
	$col = 30 if $tier =~ /(?:^|,|\b)(?:\(|--)/;
	return hl($col, $tier);
}

#######################################################################################################################################
sub extract_results_catypes {
	my @treatment_match_result = @_;
	my @results_catypes = grep { /^catype:/ } @treatment_match_result ;
	my @retval;
	for my $catype_line (sort @results_catypes ) {
		my ($fact, @tags) = Facts::string_to_fact_and_tags( $catype_line );
		push @retval, $fact;
	}
	return @retval;
}

#######################################################################################################################################
sub extract_results_prior_therapies {
	my @treatment_match_result = @_;
	my @results_prior_therapiess = grep { /^prior_therapy:/ } @treatment_match_result ;
	my @retval;
	for my $prior_therapies_line (sort @results_prior_therapiess ) {
		my ($fact, @tags) = Facts::string_to_fact_and_tags( $prior_therapies_line );
		push @retval, $fact;
	}
	return @retval;
}

#######################################################################################################################################
sub extract_results_biomarkers {
	my @treatment_match_result = @_;
	my @results_biomarker = grep { /^has_biomarker:/ } @treatment_match_result ;
	s/\s*\[.*// for @results_biomarker ;
	my @retval;
	for my $biomarker_line (sort @results_biomarker) {
		my %AMP_LOE;
		my $AMP_LOE_max = 'U';
		my %is_original_fact;
		my $biomarker = $biomarker_line;
		$biomarker =~ s/has_biomarker://;
		my %a;
		my @matched = map { s/^\Q$biomarker\E://; $_ } grep { /^\Q$biomarker\E:/ } @treatment_match_result ;
		@matched = grep { ! /^[VADSFT]$/ and ! /(?:treatment(?:_class)?:|targeted_therapy)/ and ! /pkmb_tier|sensitive_to|resistant_to/  } @matched;  
		s/^V:// for @matched ;
		@matched = uniq(@matched);
		
		my @aberration ;
		my @oncogenicity ;
		my @consequence ;
		my @specifics ;
		for my $fact_str (@matched) {
			my ($fact, @tags) = Facts::string_to_fact_and_tags( $fact_str );
			$is_original_fact{$fact} = ( exists $facts{"$biomarker:$fact"} ? 1 : 0 ) ;
			if ( $fact =~ /^\s*(?:alteration|mutation|deletion|amplification|.*expression|high|low|fusion|methylation|.*type|.*duplication)\s*$/i ) { 
				push @aberration, $fact; 
			} elsif ( $fact_str =~ /-of-function_mutation|oncogenic_mutation/ ) {
				push @consequence, $fact_str;
			} else {
				push @specifics, $fact;
			}
			if ( my @a = grep { /:oncogenicity/ } @tags ) {
				push @oncogenicity, @a;
			} 			
			my @amp_loe = sort grep { /^AMP\.LOE:/ } @tags;
			if ( scalar(@amp_loe) and $amp_loe[0] =~ /^AMP\.LOE:(.+)/ ) {
				$AMP_LOE{$fact} = $1; 
				$AMP_LOE_max = $AMP_LOE{$fact} if ( $AMP_LOE_max eq 'U' ) or (defined $AMP_LOE{$fact} and ( $AMP_LOE{$fact} lt $AMP_LOE_max ) );
# 				print STDERR "$AMP_LOE_max\t$fact\t$AMP_LOE{$fact}\n";
			} else {
				$AMP_LOE{$fact} = undef;
			}
		}
		
		@aberration = sort { ( $is_original_fact{$b} <=> $is_original_fact{$a} ) || ( ($AMP_LOE{$a} // 3) cmp ($AMP_LOE{$b} // 3) ) || ($a cmp $b ) } @aberration ;
		s/\s*\[.+?\]\s*$//   for @matched ;
		s/:oncogenicity//    for @oncogenicity ;
		s/\s*\[.+?\]\s*$//   for @oncogenicity ;
		s/:mutation_effect// for @consequence ;
		s/\s*\[.+?\]\s*$//   for @consequence ;
		@matched = grep { !/_to:/ } @matched ;
		@consequence = sort @consequence ;
		@specifics = uniq(@specifics) ;
		@oncogenicity = uniq(@oncogenicity);
		@specifics = sort { ( $is_original_fact{$b} <=> $is_original_fact{$a} ) || ( ($AMP_LOE{$a} // 3) cmp ($AMP_LOE{$b} // 3) ) || ($a cmp $b ) }  @specifics ;
		s/_/ /g for @matched;
		
# 		print STDERR "AMP_LOE_max($biomarker) = $AMP_LOE_max\n";
		my %row = (
			'biomarker'    => $biomarker ,
			'aberration'   => \@aberration,
			'specifics'    => \@specifics,
			'consequence'  => \@consequence,
			'oncogenicity' => \@oncogenicity,
			'AMP.LOE'          => \%AMP_LOE,
			'AMP.LOE.max'      => $AMP_LOE_max,
			'is_original_fact' => \%is_original_fact,
		);
		push @retval, \%row;
	}	
	@retval = sort { ( ( $$a{'AMP.LOE.max'} // '5' ) cmp ( $$b{'AMP.LOE.max'} // '5' ) ) || ($a cmp $b) } @retval;
# 	print STDERR Dumper(@retval);

	return @retval;
};

#######################################################################################################################################

sub extract_drug_sensitivity {

	my @results_to_filter = @_;
	
	@results_to_filter = sort { ( extract_pref_score($a) <=> extract_pref_score($b) ) || score_LOM($a) <=> score_LOM($b) || ($a cmp $b) } @results_to_filter;
	
	my @retval;
	for my $result_line ( @results_to_filter ) {
		my %row;
		my $ent = $result_line ;
		my $is_drug_class = ( $result_line =~ /^recommendation_tier_drug_class/ );
		$ent =~ s/^recommendation_tier\S*?://;
		$ent =~ s/\S+_(?:to|evidence)_\S+://;
		$ent =~ s/\s*\[(.*)\]\s*$//;
		my @tags = split /\s*;\s*/, $1 if defined $1 ;
		my @biomarkers = sort ( uniq( map { s/referred_from://; split /\s*\+\s*/, $_} grep { /referred_from:/ } @tags ) );
		my $biomarkers = join ", ", @biomarkers;
		my @lom = sort map { my $a = $_; $a =~ s/^LOM:\s*(\d+).*/$1/; $a } grep { /^LOM:/ } @tags; 
		my $lom = join(", ", @lom);

		my @alterations;
		
		for my $bm (@biomarkers) {
			my @matched = grep { /^\Q$bm\E:/ and ! /:treatment(?:_class)?:/ and ! /:[AVDTS](?:$|:)/ } @tags;
			if ( scalar @biomarkers == 1 ) { 
				s/^\Q$bm\E:// for @matched ;
			}
			push @alterations, @matched;
		}
		$ent =~ s/_/ /g;
		$biomarkers =~ s/_/ /g;
	
		my ($therapy, $tier) = ($ent =~ /(.*):(\S+)/);
	
		my $alterations = join ", ", @alterations ;
		$row{'alterations'} = $alterations ;
		$row{'biomarker'} = $biomarkers ;
		$row{'therapy'} = $therapy ;
		$row{'tier'} = $tier ;
		$row{'lom'} = $lom ;
		
		if ( $f_clinical_report ) {
			next if $tier eq 'R2B';
			if ( ! $is_drug_class ) {
				next if $tier =~ /^[34]/;
			}
		}
		
		push @retval, \%row;

	}
	return @retval;
}

#######################################################################################################################################

sub get_prefential_trial_id { my ($id) = ( $_[0] =~ /preferential_trial_id:(\S+)/ ); return $id; }
sub get_prefential_trial_score { my ($score) = ( $_[0] =~ /pref_trial_score:(-?[\d\.]+)(?:\b|$)/ ); return $score; }
sub cmp_preferential_trials($$) {
	my ($a, $b) = @_ ;
	return 
		( get_prefential_trial_score($a) <=> get_prefential_trial_score($b) ) ||
		( get_prefential_trial_id($a) cmp get_prefential_trial_id($b) );
}

sub extract_preferential_trials {
	my @results_preferential_trials = @_;
	my @retval;
	for my $line ( @results_preferential_trials ) {
		my ($trial_id, $tags) = ( $line =~ /^preferential_trial_id:(\S+)\s*\[(.+)\]\s*$/ );
		my %tags = Facts::string_to_tags_hashed($line);
		my ($fact, @tags) = Facts::string_to_fact_and_tags($line);
		
		my %referring_catypes = map { s/catype://g; $_ => 1 } ( grep { /^\(?\s*catype:(.+)/ } @tags );
		my %referring_genes = map { $_ => 1 } ( $tags =~ /referred_from:(.+?)[;\]]/g );
		
		my %referring_drugs        = map { /recommendation_tier:(.+):\S+/; $1 => 1 }            grep { /recommendation_tier:.+:\S+/ }            @tags ;
		my %referring_drug_classes = map { /recommendation_tier_drug_class:(.+):\S+/; $1 => 1 } grep { /recommendation_tier_drug_class:.+:\S+/ } @tags ;
		
		my %referring_alterations; 
		
		for my $biomarker ( keys %referring_genes ) {
			$referring_alterations{$_}++ for grep { /^$biomarker:(.+)/ } @tags;
		}
		
		my $full_title         = Facts::unescape( $tags{'studytitle'}        // '' );
		my $trialacronym       = Facts::unescape( $tags{'trialacronym'}      // '' );
		my $trialphase         = Facts::unescape( $tags{'phase'}      // '' );
		my $recruitmentstatus  = Facts::unescape( $tags{'recruitmentstatus'} // '' );
		my $matched_drug_names = Facts::unescape( $tags{'drug_list'}         // '' );
		my $matched_trial_tier                 =  $tags{'transitive_efficacy_tier'};
		my $matched_trial_class_tier           =  $tags{'transitive_class_efficacy_tier'};
		my $matched_trial_drug_maturity        =  $tags{'drug_maturity_tier'};            
		my $matched_trial_drug_class_maturity  =  $tags{'drug_class_maturity_tier'}; 
		my $matched_trial_combo_maturity       =  $tags{'combo_maturity_tier'};   
		my $matched_trial_combo_class_maturity =  $tags{'combo_class_maturity_tier'};
		my $matched_trial_phase_tier           =  $tags{'trial_phase_tier'};
		my $matched_trial_biomarker_tier       =  $tags{'biomarker_tier'};
		my $matched_trial_match_criteria_score =  $tags{'trial_match_criteria_score'};
		my $matched_trial_referring_drug_classes_score =  $tags{'referring_drug_classes_score'};
		my $matched_trial_LOM_score            =  $tags{'LOM'} // 'Not assessed';
		my $ext_weblink                        =  $tags{'ext_weblink'};
		my $notes                              =  join("; ", grep { /^\*/ } @tags);
		my $preferential_trial_score           =  $tags{'pref_trial_score'}  // '';
		
		my @matched_drug_names   = split /\s*;\s*/, $matched_drug_names ;
		@matched_drug_names = uniq( map { Therapy::get_preferred_drug_name($_) } @matched_drug_names );
		@matched_drug_names = sort @matched_drug_names ;
		
		my $matched_drug_classes = Facts::unescape( $tags{'drug_classes'} );
		my @matched_drug_classes = split /\s*;\s*/, $matched_drug_classes ;
# 		@matched_drug_classes    = map { exists $referring_drug_classes{$_} ? "<b>$_</b>" : $_  } grep { ! / \+ / } @matched_drug_classes ;
		my @healthconditions   = split /;\s*/, ( Facts::unescape( $tags{'healthcondition'} // '' ) ); 
		my @postcodes          = sort map { s/^\s*-\s*//; $_} split(/;\s*/, ( Facts::unescape( $tags{'postcode'} // '' ) )); 
		
		my %trial_match_criteria = map { $_ => 1 } ( $tags =~ /trial_match_criteria:(.+?)[;\]]/g );
		my @trial_match_criteria = sort map { s/_/ /g; $_ } keys %trial_match_criteria ;
		
		my %row = (
			'trial_id' => $trial_id,
			'full_title' => $full_title,
			'trialacronym' => $trialacronym,
			'trialphase' => $trialphase,
			'recruitmentstatus' => $recruitmentstatus,
			'referring_catypes' => \%referring_catypes,
			'referring_genes' => \%referring_genes,
			'referring_alterations' => \%referring_alterations,
			'referring_drugs' => \%referring_drugs,
			'referring_drug_classes' => \%referring_drug_classes,
			'matched_drug_names' => \@matched_drug_names,
			'matched_drug_classes' => \@matched_drug_classes,
			'matched_trial_tier' => $matched_trial_tier,
			'matched_trial_class_tier' => $matched_trial_class_tier,
			'matched_trial_biomarker_tier' => $matched_trial_biomarker_tier,
			'matched_trial_drug_maturity' => $matched_trial_drug_maturity,
			'matched_trial_drug_class_maturity' => $matched_trial_drug_class_maturity,
			'matched_trial_combo_maturity' => $matched_trial_combo_maturity,
			'matched_trial_combo_class_maturity' => $matched_trial_combo_class_maturity,
			'matched_trial_match_criteria_score' => $matched_trial_match_criteria_score,
			'matched_trial_referring_drug_classes_score' => $matched_trial_referring_drug_classes_score,
			'matched_trial_phase_tier_score' => $matched_trial_phase_tier,
			'trial_match_criteria' => \@trial_match_criteria,
			'healthcondition' => \@healthconditions,
			'postcodes' => \@postcodes,
			'LOM' => $matched_trial_LOM_score,
			'ext_weblink' => $ext_weblink,
			'notes' => $notes,
			'pref_trial_score' => $preferential_trial_score
		);
		push @retval, \%row;
	}
	return @retval;
}

#######################################################################################################################################
sub gen_catypes_terminal {
	my $output = hl(37, "CANCER TYPES MATCHED")."\n";
	for my $row ( extract_results_catypes(@treatment_match_result) ) {
		$row =~ s/catype.*?://;
		$output .= "$row\n";
	}
	return $output;
}

sub gen_catypes_HTML {
	my $output = "<h2>CANCER TYPES MATCHED</h2>\n";
	for my $row ( extract_results_catypes(@treatment_match_result) ) {
		$row =~ s/catype.*?://;
		$output .= "$row<br>\n";
	}
	return $output;
}

#######################################################################################################################################
sub gen_prior_therapies_terminal {
	my $output = hl(37, "PRIOR THERAPIES MATCHED")."\n";
	for my $row ( extract_results_prior_therapies(@treatment_match_result) ) {
		$row =~ s/prior_therap.+?://;
		$output .= "$row\n";
	}
	return $output;
}

sub gen_prior_therapies_HTML {
	my $output = "<h2>PRIOR THERAPIES MATCHED</h2>\n";
	for my $row ( extract_results_prior_therapies(@treatment_match_result) ) {
		$row =~ s/prior_therap.+?://;
		$output .= "$row<br>\n";
	}
	return $output;
}


#######################################################################################################################################

sub gen_biomarker_terminal_append_AMP_LOE {
	my $x = shift;
	my $bm_hashref = shift;
	return undef if ! defined $x;
	my $retval = ( exists $$bm_hashref{'is_original_fact'}{$x} and defined $$bm_hashref{'is_original_fact'}{$x} and length $$bm_hashref{'is_original_fact'}{$x} and $$bm_hashref{'is_original_fact'}{$x} ) ? 
		hl(37, $x) : $x;
# 	print ">>>> $x\t$$bm_hashref{'AMP.LOE'}{$x}\n";
	$retval .= " [AMP:".thl_AMP($$bm_hashref{'AMP.LOE'}{$x})."]" if ( exists $$bm_hashref{'AMP.LOE'}{$x} and defined $$bm_hashref{'AMP.LOE'}{$x} and length $$bm_hashref{'AMP.LOE'}{$x} );
	return $retval;
}

sub gen_biomarker_terminal {
	my $output = hl(37, "BIOMARKER ALTERATIONS")."\n";
	
	my @lines ;
	push @lines, join("\t", "Gene", "Aberration", "Specification", "Consequence", "Oncogenicity")."\n";
	for my $row ( extract_results_biomarkers(@treatment_match_result) ) {
		for my $i (0 .. 20) {
			my $line  = '';
			$line .= ( ($i == 0) ? $$row{'biomarker'} : '' )."\t";
			$line .= ( gen_biomarker_terminal_append_AMP_LOE( @{ $$row{'aberration'} }[$i]    , $row)  // '' )."\t";
			$line .= ( gen_biomarker_terminal_append_AMP_LOE( @{ $$row{'specifics'} }[$i]     , $row)  // '' )."\t";
			$line .= ( gen_biomarker_terminal_append_AMP_LOE( @{ $$row{'consequence'} }[$i]   , $row)  // '' )."\t";
			$line .= ( gen_biomarker_terminal_append_AMP_LOE( @{ $$row{'oncogenicity'} }[$i]  , $row)  // '' );
			$line .= "\n";
			last if $line =~ /^\s*$/;
			push @lines, $line ;
		}
	}
	$output .= $_ for column(@lines) ;
	return $output ;
}

sub gen_biomarker_report {
	my $output ; # = "List of alterated genes, biomarkers, or tumoural molecular features\n";
	
	$output .= join("\t", "Gene", "Aberration", "Specification", "Consequence", "Oncogenicity")."\n";
	for my $row ( extract_results_biomarkers(@treatment_match_result) ) {
		$output .= join("\t", $$row{'biomarker'} , 
			join(";", @{ $$row{'aberration'} } ) ,
			join(";", @{ $$row{'specifics'} } ) ,
			join(";", @{ $$row{'consequence'} } ) ,
			join(";", @{ $$row{'oncogenicity'} } ) ,
			)."\n";
	}
	return $output ;
}

sub clean_HTML {
	my $x = shift;
	$x =~ s/&/&amp;/g;
	$x =~ s/>/&gt;/g;
	$x =~ s/</&lt;/g;
	return $x;
}

sub gen_HTML_biomarker_report {
	my $output = "<h2>List of alterated genes, biomarkers, or tumoural molecular features</h2>\n<table id=trial_list>\n";
	$output .= "<tr>".join("", ( map { "<th>$_</th>" } ("Gene", "Aberration", "Specification", "Consequence", "Oncogenicity") ) )."</tr>\n"; # "Status", "Relevance", "Hospital", 
	my $cnt = 0;
	for my $row ( extract_results_biomarkers(@treatment_match_result) ) {
		my $oe = (($cnt % 2) ? 'e': 'o');
		$output .= "<tr class=$oe>".join("", (map {"<td>$_</td>"} (
			$$row{'biomarker'} , 
			join("<br/>\n", @{ $$row{'aberration'} } ) ,
			join("<br/>\n", ( map { clean_HTML($_) } @{ $$row{'specifics'} } ) ) ,
			join("<br/>\n", @{ $$row{'consequence'} } ) ,
			join("<br/>\n", @{ $$row{'oncogenicity'} } ) ,
			)
		) )."</tr>\n";
		++$cnt ;
	}
	
	$output .= "</table>\n";
	$output .= "$cnt entities defined.<br>";
	return $output ;
}

#######################################################################

sub gen_drug_sens_terminal {
	my $output = "\e[1;37mTHERAPY SENSITIVITY/RESISTANCE\e[0m\n";
	
	our @title_rep ; 
	
	sub gen_drug_sens_table_terminal {
		my $title = shift;
		my @results_to_filter = @_;
		my @lines ;
		push @title_rep, $title;
		push @lines, "___$#title_rep\n";
		push @lines, join("\t", "Therapy/Class", 'Tier', 'Biomarker(s)', 'LOM', 'Alteration(s)')."\n";
		for my $row ( extract_drug_sensitivity(@results_to_filter) ) {
			my $f_hl = 0;
			$f_hl = 37 if $$row{'tier'} =~ /^1|^2/;
			$f_hl = 30 if $$row{'tier'} =~ /^3B|^4|U|^R2B/;
			$f_hl = 31 if $$row{'tier'} =~ /^R1/;
			push @lines, join("\t", ($f_hl ? hl($f_hl, $$row{'therapy'}) : $$row{'therapy'}) , thl( $$row{'tier'} ), $$row{'biomarker'}, $$row{'lom'}, $$row{'alterations'})."\n";
		}
		return @lines ;
	}
	my @lines;
	push @lines, gen_drug_sens_table_terminal("Sensitive to therapies or therapy combinations", @results_sens_to_drugs )  if scalar @results_sens_to_drugs ;
	push @lines, gen_drug_sens_table_terminal("Sensitive to therapy classes", @results_sens_drug_class)                   if scalar @results_sens_drug_class ;
	push @lines, gen_drug_sens_table_terminal("Resistant to therapies or therapy combinations", @results_resi_to_drugs )  if scalar @results_resi_to_drugs ;
	push @lines, gen_drug_sens_table_terminal("Resistant to therapy classes", @results_resi_drug_class)                   if scalar @results_resi_drug_class ;
	@lines = column(@lines) ;
	@lines = map { s/___(\d+)/\n$title_rep[$1]/g; $_ } grep { defined } @lines ;
	$output .= $_ for @lines ;
	return $output ;
}


sub gen_drug_sens_report {
	my $output ; # = "List of sensitived and resistant therapy and therapy classes\n";

	sub gen_drug_sens_table {
		my $title = shift;
		my @results_to_filter = @_;
		my $output = '';
		$output .= "$title\n".
			join("\t", 'Therapy or class', 'Tier', 'Biomarker(s)', 'LOM', 'Alteration(s)')."\n";
		
		for my $row ( extract_drug_sensitivity(@results_to_filter) ) {
			$output .= join("\t", $$row{'therapy'}, $$row{'tier'}, $$row{'biomarker'}, $$row{'lom'}, $$row{'alterations'})."\n";
		}
		return $output ;
	}

	$output .= gen_drug_sens_table("Sensitive to therapies or therapy combinations", @results_sens_to_drugs );
	$output .= gen_drug_sens_table("Sensitive to therapy classes", @results_sens_drug_class);
	$output .= gen_drug_sens_table("Resistant to therapies or therapy combinations", @results_resi_to_drugs );
	$output .= gen_drug_sens_table("Resistant to therapy classes", @results_resi_drug_class);

	return $output;
}


sub gen_HTML_drug_sens_report {
	my $output = "<h2>List of sensitive and resistant therapy and therapy classes</h2>\n";

	sub gen_HTML_drug_sens_table {
		my $title = shift;
		my @results_to_filter = @_;
		my $output = '';
		$output .= "<h3>$title</h3>";
		$output .= "<table style='border-collapse:collapse' id=trial_list>\n";
		$output .= "<tr>".join("", ( map { "<th>$_</th>" } ('Therapy or class', 'Tier', 'Biomarker(s)', 'LOM', 'Alteration(s)') ) )."</tr>\n"; # "Status", "Relevance", "Hospital", 
		my $cnt = 0;
		
		for my $row ( extract_drug_sensitivity(@results_to_filter) ) {
			my $oe = (($cnt % 2) ? 'e': 'o');
			
			my $bgcolor = '';
			$bgcolor = '#3deb3d' if $$row{'tier'} =~ /^1/;
			$bgcolor = '#e6ff00' if $$row{'tier'}=~ /^1B|1R/;
			$bgcolor = '#eeee99' if $$row{'tier'}=~ /^2/;
			$bgcolor = '#aaccff' if $$row{'tier'}=~ /^3/;
			$bgcolor = '#eeeeee' if $$row{'tier'}=~ /^3B/;
			$bgcolor = '#cccccc' if $$row{'tier'}=~ /^4/;
			$bgcolor = '#ff9999' if $$row{'tier'}=~ /^R1/;
			$bgcolor = '#ffcc99' if $$row{'tier'}=~ /^R2/;
			
			$output .= "<tr class=$oe style='background:$bgcolor'>".join("", (map {"<td style='border:1px solid #888888'>$_</td>"} ( "<b>$$row{'therapy'}</b>", $$row{'tier'}, $$row{'biomarker'}, $$row{'lom'}, clean_HTML($$row{'alterations'}) ) ) )."</tr>\n";
			++$cnt ;
		}
		$output .= "</table>"; # ."</td>\n";
		return $output ;
	}

	$output .= "<table style='border:none'><tr>\n";
	$output .= "<td style='vertical-align:top;'>\n";
	$output .= gen_HTML_drug_sens_table("Sensitive to therapies or therapy combinations", @results_sens_to_drugs );
	$output .= "<td/>\n";
	$output .= "<td style='vertical-align:top;'>\n";
	$output .= gen_HTML_drug_sens_table("Sensitive to therapy classes", @results_sens_drug_class);
	$output .= "<td/>\n";
	$output .= "<tr></tr>\n";
	$output .= "<td style='vertical-align:top;'>\n";
	$output .= gen_HTML_drug_sens_table("Resistant to therapies or therapy combinations", @results_resi_to_drugs );
	$output .= "<td/>\n";
	$output .= "<td style='vertical-align:top;'>\n";
	$output .= gen_HTML_drug_sens_table("Resistant to therapy classes", @results_resi_drug_class);
	$output .= "<td/>\n";
	$output .= "</tr></table>\n";

	return $output;
}

sub mkhref_evidence {
	my $x = shift;
	for ($x) {
		/^(\d+)$/    and do { return "<a href='https://pubmed.ncbi.nlm.nih.gov/$1' target=_blank>$1</a>"; last; };
		/^(10\..+)$/ and do { return "<a href='https://doi.org/$1' target=_blank>$1</a>"; last; };
		/^(NCT\d+)$/ and do { return "<a href='https://clinicaltrials.gov/ct2/show/NCT03798626' target=_blank>"; last; };
	}
	return $x;
}

sub sort_results_treatment {
	# Sort by curated tier first, then sort by level of matching, then alphabetically
	@results_treatment = sort { ( score_tier($a) <=> score_tier($b) ) || ( score_LOM($a) <=> score_LOM($b) ) || ($a cmp $b) } @results_treatment;
}

sub gen_HTML_biomarker_evidence_report {
	my $output = "<h2>Evidence of matched target therapies:</h2>\n<table id=trial_list>\n";
	$output .= "<tr>".join("", ( map { "<th>$_</th>" } ("Biomarker", "Drug", "Level of Evidence", "Reference") ) )."</tr>\n"; # "Status", "Relevance", "Hospital", 

	sort_results_treatment();

	my $cnt = 0;
	for my $treatment_line (@results_treatment) {
		my ($kb_loe)     = ( $treatment_line =~ /\(\w+ LOE: (.+?)\)/ );
		my %tags = Facts::string_to_tags_hashed( $treatment_line );
		my @evidence = ( exists $tags{'evidence'} ? ( split /\s*[;,]\s*/, $tags{'evidence'}) : () );
		@evidence = map { mkhref_evidence($_) } @evidence ;
		$treatment_line =~ s/(?:\(\w+|\[LOM).*//;
		my ($gene_trigger, $dummy, $treatment)  =  split /:/, $treatment_line ;
		my $oe = (($cnt % 2) ? 'e': 'o');
		$kb_loe .= " [repurposed]" if $kb_loe =~ /, in/;
# 		next if $kb_loe =~ /, in/;
		
		$output .= "<tr class=$oe>".join("", (map {"<td>$_</td>"} (
			$gene_trigger, 
			$treatment, 
			$kb_loe,
			join(', ', @evidence)
			)
		) )."</tr>\n";
	}
	$output .= "</table><br/>\n";

# 	print $output ;
	return $output;
}

sub gen_biomarker_evidence_report {
	my $output ; # = "Evidence of matched target therapies:\n";
	$output .= join("\t", "Biomarker", "Drug", "Level of Evidence", "Reference")."\n"; # "Status", "Relevance", "Hospital", 

	sort_results_treatment();

	my $cnt = 0;
	for my $treatment_line (@results_treatment) {
		my ($kb_loe)     = ( $treatment_line =~ /\(\w+ LOE: (.+?)\)/ );
		my %tags = Facts::string_to_tags_hashed( $treatment_line );
		my @evidence = ( exists $tags{'evidence'} ? ( split /\s*[;,]\s*/, $tags{'evidence'}) : () );
		$treatment_line =~ s/(?:\(\w+|\[LOM).*//;
		my ($gene_trigger, $dummy, $treatment)  =  split /:/, $treatment_line ;
# 		next if $kb_loe =~ /, in/;
		$kb_loe .= " [repurposed]" if $kb_loe =~ /, in/;
		
		$output .= join("\t", $gene_trigger,  $treatment,  $kb_loe, join(', ', @evidence) )."\n";
	}
	$output .= "\n";
	return $output;
}

sub gen_biomarker_evidence_terminal {
	my $output = "\e[1;37mEVIDENCE OF BIOMARKER-MATCHED THERAPIES\e[0m\n";
	my @lines ;
	push @lines, join("\t", "Biomarker", "Drug", "KB Name", "Level of Evidence", "Reference")."\n"; # "Status", "Relevance", "Hospital", 

	sort_results_treatment();

	my $cnt = 0;
	for my $treatment_line (@results_treatment) {
		my ($kb_name, $kb_loe)     = ( $treatment_line =~ /\((\S+) LOE: (.+?)\)/ );
		my %tags = Facts::string_to_tags_hashed( $treatment_line );
		my @evidence = ( exists $tags{'evidence'} ? ( split /\s*[;,]\s*/, $tags{'evidence'}) : () );
		$treatment_line =~ s/(?:\(\w+|\[LOM).*//;
		my ($gene_trigger, $dummy, $treatment)  =  split /:/, $treatment_line ;
# 		next if $kb_loe =~ /, in/;
		$kb_loe .= " [repurposed]" if $kb_loe =~ /, in/;
		
		push @lines, join("\t", $gene_trigger,  $treatment,  $kb_name, thl( $kb_loe ), join(', ', @evidence) )."\n";
	}
	@lines = column(@lines) ;
	$output .= join '', @lines;
	return $output;
}

#############################################################################################################################################################

sub extract_therapy_summaries {
	my @retval;

	for my $line ( @results_therapy_recommendations ) {
		my ($therapy_recommendation_str, @tags) = Facts::string_to_fact_and_tags( $line );
		next if ( $therapy_recommendation_str !~ /^therapy_recommendation:(.+):([[:alnum:]]+)$/ );
		my ($therapy, $tier) = ($1, $2);
		
		my %tags = Facts::string_to_tags_hashed( $line );
		
		my @evidence = map { my $a = $_; $a =~ s/^evidence:\s*//; my @a = split /\s*[,;]\s*/, $a; @a } grep { /^evidence:/ } @tags; 
		my @referred_from = sort map { my $a = $_; $a =~ s/^referred_from:\s*//; my @a = split /\s*[,;]\s*/, $a; @a } grep { /^referred_from:/ } @tags; 
		my @lom = sort map { my $a = $_; $a =~ s/^LOM:\s*(\d+).*/$1/; $a } grep { /^LOM:/ } @tags; 

		my %row = (
			'therapy'         => $therapy,
			'tier'            => $tier,
			'tier_drug_class' => $tags{'therapy_recommendation_tier_drug_class'},
			'rec_score'       => $tags{'therapy_recommendation_score'},
			'referred_from'   => join(", ", @referred_from),
			'evidence'        => join(", ", @evidence),
			'lom'             => join(", ", @lom)
		);
		
		push @retval, \%row;
	}
	return @retval;
}


sub gen_therapy_summary_report {
	my $output = join("\t", "Therapy", "Tier", "Tier (drug class)", "Biomarker", "LOM", "Evidence")."\n";
	my @lines;

	for my $row ( sort { $$a{'rec_score'} <=> $$b{'rec_score'} } extract_therapy_summaries(@results_therapy_recommendations) ) {
		push @lines, join("\t", $$row{'therapy'}, $$row{'tier'}, $$row{'tier_drug_class'}, $$row{'referred_from'}, $$row{'lom'}, $$row{'evidence'});
	}

	$output .= join("\n", @lines)."\n";
	return $output;
}

sub gen_therapy_summary_terminal {
	my $output = "\e[1;37mTHERAPY RECOMMENDATIONS\e[0m\n";
	
	my @lines;
	
	push @lines, join("\t", "Therapy", "Tier", "Tier (drug class)", "Biomarker", "LOM", "Evidence")."\n";
	
	for my $row ( sort { $$a{'rec_score'} <=> $$b{'rec_score'} } extract_therapy_summaries(@results_therapy_recommendations) ) {
		my $f_hl;
		$f_hl = 37 if $$row{'tier'} =~ /^1|^2/;
		$f_hl = 30 if $$row{'tier'} =~ /^3B|^4|^R2B/;
		$f_hl = 31 if $$row{'tier'} =~ /^R1/;
# 		push @lines, join( "\t", ( map { $$row{$_} }  qw(therapy tier tier_drug_class evidence) ) )."\n";
		push @lines, join("\t", ($f_hl ? hl($f_hl, $$row{'therapy'}) : $$row{'therapy'}) , thl( $$row{'tier'} ), thl( $$row{'tier_drug_class'} ), $$row{'referred_from'}, $$row{'lom'}, $$row{'evidence'})."\n";
	}
	
	@lines = column(@lines) ;
	$output .= $_ for @lines ;
	return $output ;
}

sub gen_HTML_therapy_summary_report {
	my $output = "<h2>Therapy recommendations</h2>\n";

	my $title = shift;
	$output .= "<table style='border-collapse:collapse' class=trial_list>\n";
	$output .= "<tr>".join("", ( map { "<th>$_</th>" } ("Therapy", "Tier", "Tier (drug class)", "Biomarker", "LOM", "Evidence") ) )."</tr>\n"; # "Status", "Relevance", "Hospital", 
	my $cnt = 0;
	
	for my $row ( sort { $$a{'rec_score'} <=> $$b{'rec_score'} } extract_therapy_summaries(@results_therapy_recommendations) ) {
		my $oe = (($cnt % 2) ? 'e': 'o');
		
		my $bgcolor = '';
		$bgcolor = '#3deb3d' if $$row{'tier'} =~ /^1/;
		$bgcolor = '#e6ff00' if $$row{'tier'}=~ /^1B|1R/;
		$bgcolor = '#eeee99' if $$row{'tier'}=~ /^2/;
		$bgcolor = '#aaccff' if $$row{'tier'}=~ /^3/;
		$bgcolor = '#eeeeee' if $$row{'tier'}=~ /^3B/;
		$bgcolor = '#cccccc' if $$row{'tier'}=~ /^4|U/;
		$bgcolor = '#ff9999' if $$row{'tier'}=~ /^R1/;
		$bgcolor = '#ffcc99' if $$row{'tier'}=~ /^R2/;
		
		my @evidence = ( exists $$row{'evidence'} ? ( split /\s*[;,]\s*/, $$row{'evidence'}) : () );
		@evidence = sort map { mkhref_evidence($_) } @evidence ;
		
		$output .= "<tr class=$oe style='background:$bgcolor'>".join("", (
			map {"<td style='border:1px solid #888888'>$_</td>"} ( 
			"<b>$$row{'therapy'}</b>", 
			$$row{'tier'}, 
			$$row{'tier_drug_class'}, 
			$$row{'referred_from'}, 
			$$row{'lom'}, 
			join(", ", @evidence) ) )
			)."</tr>\n";
		++$cnt ;
	}
	$output .= "</table>"; # ."</td>\n";

	return $output;
}

#############################################################################################################################################################

sub gen_preferential_trial_report {
	@results_preferential_trials = sort { cmp_preferential_trials($a, $b) } @results_preferential_trials ;

	my $output ; # = "List of biomarker matched clinial trials:\n";
	$output .= join("\t", "Rank", "Trial ID", "Preferential Trial Score", 
		"Transitive Class Efficacy", "Transitive Efficacy", "Drug maturity", "Drug class maturity", "Combo maturity", "Combo class maturity", "Biomarker tier score", "Trial phase tier score", "Trial match criteria", "Level of matching",
		"Drugs", "Drug classes", "Cancer types", "Phase", "Full title", "Postcode", "External weblink", "Assumptions"
		)."\n";

	my $cnt = 1;

	for my $row ( extract_preferential_trials(@results_preferential_trials) ) {
		my @matched_drug_names = @{ $$row{'matched_drug_names'} };
		@matched_drug_names = @matched_drug_names ;
		my $matched_drug_names = join("; ", @matched_drug_names);

		my @matched_drug_classes = @{ $$row{'matched_drug_classes'} };
		@matched_drug_classes    = @matched_drug_classes ;
		my $matched_drug_classes = join("; ", @matched_drug_classes);
		
		my $healthcondition       = join("; ", @{ $$row{'healthcondition'} } ); 
		my $postcodes             = join("; ", @{ $$row{'postcodes'} });
		my $trial_match_criteria  = join("; ", @{ $$row{'trial_match_criteria'} });

		$output .= join("\t", (
				$cnt, 
				$$row{'trial_id'}, 
				$$row{'pref_trial_score'},
				$$row{'matched_trial_class_tier'},
				$$row{'matched_trial_tier'},
				$$row{'matched_trial_drug_maturity'},
				$$row{'matched_trial_drug_class_maturity'},
				$$row{'matched_trial_combo_maturity'},
				$$row{'matched_trial_combo_class_maturity'},
				$$row{'matched_trial_phase_tier_score'},
				$$row{'matched_trial_biomarker_tier'},
				$trial_match_criteria,
				$$row{'LOM'},
				$matched_drug_names,
				$matched_drug_classes,
				$healthcondition,
				$$row{'trialphase'},				
				(length( $$row{'$trialacronym'} ) ? "$$row{'$trialacronym'} - " : "").$$row{'full_title'},
				$postcodes,
				$$row{'ext_weblink'},
				$$row{'notes'}
				)
			)."\n";
		++$cnt ;
	}

	return $output ;
}


sub gen_preferential_trial_terminal {
	@results_preferential_trials = sort { cmp_preferential_trials($a, $b) } @results_preferential_trials ;

	my $output = "\n\n\e[1;37mBIOMARKER-MATCHED, RATIONALLY PRIORITISED CLINICAL TRIALS\e[0m\n";
	
	my @lines ;

	my $cnt = 1;

	for my $row ( extract_preferential_trials(@results_preferential_trials) ) {
		my @matched_drug_names = grep { ! / \+ / } @{ $$row{'matched_drug_names'} };
		@matched_drug_names = @matched_drug_names ;
		my $matched_drug_names = join("; ", @matched_drug_names);

		my @matched_drug_classes = grep { ! / \+ / } @{ $$row{'matched_drug_classes'} };
		@matched_drug_classes    = @matched_drug_classes ;
		my $matched_drug_classes = join("; ", @matched_drug_classes);
		
		my $healthcondition      = join("; ", @{ $$row{'healthcondition'} } ); 
		my $postcodes            = join("; ", @{ $$row{'postcodes'} });
		my $trial_match_criteria = join("; ", @{ $$row{'trial_match_criteria'} });

		push @lines, join("\t", 'Rank',                       $cnt)."\n";
		push @lines, join("\t", 'Trial ID',                   hl(37, $$row{'trial_id'}) )."\n";
		push @lines, join("\t", "Full title",                 (length( $$row{'$trialacronym'} ) ? "$$row{'$trialacronym'} - " : "").$$row{'full_title'} )."\n";
		if (! $f_list_trials_terse ) {
			push @lines, join("\t", "Trial phase", $$row{'trialphase'} )."\n";
			push @lines, join("\t", "Drugs",                      $matched_drug_names )."\n";
			push @lines, join("\t", "Drug classes",               $matched_drug_classes )."\n";
			push @lines, join("\t", "Transitive class efficacy tier",  thl( $$row{'matched_trial_class_tier'} ) )."\n";
			push @lines, join("\t", "Transitive efficacy tier",   thl( $$row{'matched_trial_tier'} ) )."\n";
			push @lines, join("\t", "Trial phase tier",           thl( $$row{'matched_trial_phase_tier_score'} ) )."\n";
			push @lines, join("\t", "Drug maturity tier",         thl( $$row{'matched_trial_drug_maturity'} ) )."\n";
			push @lines, join("\t", "Drug class maturity tier",   thl( $$row{'matched_trial_drug_class_maturity'} ) )."\n";
			push @lines, join("\t", "Combo maturity tier",        thl( $$row{'matched_trial_combo_maturity'} ) )."\n";
			push @lines, join("\t", "Combo class maturity tier",  thl( $$row{'matched_trial_combo_class_maturity'} ) )."\n";
			push @lines, join("\t", "Referred drug classes score", $$row{'matched_trial_referring_drug_classes_score'} )."\n";
			push @lines, join("\t", "Trial match criteria score", $$row{'matched_trial_match_criteria_score'} )."\n";
			push @lines, join("\t", "Biomarker tier",             thl_AMP( $$row{'matched_trial_biomarker_tier'} ) )."\n";
			push @lines, join("\t", "Referring cancer types",     join("; ", sort( uniq( keys %{ $$row{'referring_catypes'} } ) ) ) )."\n";
			push @lines, join("\t", "Referring biomarker",        join("; ", sort( uniq( keys %{ $$row{'referring_genes'} } ) ) ) )."\n";
			push @lines, join("\t", "Referring alterations",      join("; ", sort( uniq( keys %{ $$row{'referring_alterations'} } ) ) ) )."\n";
			push @lines, join("\t", "Referring drugs",            join("; ", sort( uniq( keys %{ $$row{'referring_drugs'} } ) ) ) )."\n";
			push @lines, join("\t", "Referring drug classes",     join("; ", sort( uniq( keys %{ $$row{'referring_drug_classes'} } ) ) ) )."\n";
		}
		push @lines, join("\t", "Trial match criteria",       $trial_match_criteria )."\n";
		push @lines, join("\t", "Level of matching",          $$row{'LOM'} )."\n";
		push @lines, join("\t", "Health conditions",          $healthcondition )."\n";
		push @lines, join("\t", "Postcodes",                  $postcodes )."\n";
		push @lines, join("\t", "External weblink",           $$row{'ext_weblink'} )."\n";
		if ( length $$row{'notes'} ) {
			my $assumptions_cnt = 0;
			my @assumptions = split /\s*;\s*/, $$row{'notes'};
			for my $assumption (@assumptions) {
				push @lines, join("\t", ($assumptions_cnt ? "" : "Assumptions"), hl(31, $assumption) )."\n" ;
				++$assumptions_cnt ;
			}
		}
		push @lines, join("\t", " " )."\n";
		++$cnt ;
	}
	
	@lines = column_wrap_text (80, @lines); 
	
# 	print map {"$_\n" } @lines;
	@lines = grep { defined } column (@lines);
	$output .= join "", @lines;
	
	return $output ;
}



sub mk_trial_href {
	my $x = shift;
	
	$x =~ s!\b((?:NCT|ACTRN)\d+[Pp]?)\b!'<a href="'.ClinicalTrials::get_trial_href($1).'" target=_blank>'.$1.'</a>'!ge ;
	
	return $x;
}


sub gen_HTML_preferential_trial_report {
	@results_preferential_trials = sort { cmp_preferential_trials($a, $b) } @results_preferential_trials ;

	my $output = "<h2>List of biomarker matched clinial trials:</h2>\n";
	$output .= "<table id=trial_list>\n";
	$output .= "<tr>".join("", ( map { "<th>$_</th>" } ("Trial ID",  "Transitive Class Efficacy", "Transitive Efficacy", "Drug maturity", "Drug class maturity", "Combo maturity", "Combo class maturity", "Trial match criteria",
		"Drugs", "Drug classes", "Cancer types", "Full title",  "Postcode", "Notes") ) )."</tr>\n"; # "Status", "Relevance", "Hospital", 

	my $cnt = 0;

	for my $row ( extract_preferential_trials(@results_preferential_trials) ) {
		my @matched_drug_names = @{ $$row{'matched_drug_names'} };
		@matched_drug_names = map { exists $$row{'referring_drugs'}{$_} ? "<b>$_</b>" : $_  } grep { ! / \+ / } @matched_drug_names ;
		my $matched_drug_names = join("<br/>", @matched_drug_names);

		my @matched_drug_classes = @{ $$row{'matched_drug_classes'} };
		@matched_drug_classes    = map { exists $$row{'referring_drug_classes'}{$_} ? "<b>$_</b>" : $_  } grep { ! / \+ / } @matched_drug_classes ;
		my $matched_drug_classes    = join("<br/>", @matched_drug_classes);
		
		my $healthcondition      = join(",<br/>", @{ $$row{'healthcondition'} } ); 
		my $postcodes            = "<div style='column-count:2'>".join("<br/>\n", @{ $$row{'postcodes'} })."</div>\n";
		my $trial_match_criteria = join("; ", @{ $$row{'trial_match_criteria'} });
		my $notes                = $$row{'notes'} ;

		my $oe = (($cnt % 2) ? 'e': 'o');
		$output .= "<tr class=$oe>".join("", ( map {"<td>$_</td>"} (
				mk_trial_href($$row{'trial_id'}), 
	# 			$$row{'$recruitmentstatus'},
				$$row{'matched_trial_class_tier'},
				$$row{'matched_trial_tier'},
				$$row{'matched_trial_drug_maturity'},
				$$row{'matched_trial_drug_class_maturity'},
				$$row{'matched_trial_combo_maturity'},
				$$row{'matched_trial_combo_class_maturity'},
				$trial_match_criteria ,
				$matched_drug_names,
				$matched_drug_classes,
				$healthcondition,
				(length( $$row{'$trialacronym'} ) ? "<b>[$$row{'$trialacronym'}]</b> " : "").$$row{'full_title'},
				$postcodes,
				$notes
				)
			)
		)."</tr>\n";
		++$cnt ;
	}
	$output .= "</table>\n";
	$output .= "$cnt trials identified.<br>\n";
	return $output ;
}




if ($f_webpage) {
	$f_title //= '';
	print <<HTMLHEAD
<html> <head> <title>$POTTR_title - $f_title </title> </head>
<style>
body {font-family:roboto slab;}
textarea#variants {width:450px; height:200px; }
span.c31m { color:#880000; } 
span.c32m { color:#008800; } 
span.c33m { color:#888800; } 
span.c36m { color:#008888; }
td { max-width: 25%; } 
h1 {font-size: 24pt;}
table#trial_list      { border:1px solid black; }
table#trial_list th   { background: black; color:white; max-width: 25%; }
table#trial_list td   { max-width: 600px;  word-wrap:break-word; vertical-align: top;} 
table#trial_list tr.o { background: #ffffcc; }
table#trial_list tr.e { background: #ffffff; }
div.debug  {background:#cccccc;}
</style>
<body>
<h1>$POTTR_title - $f_title </h1>
HTMLHEAD
;
	print gen_catypes_HTML();
	print gen_prior_therapies_HTML();
	print gen_HTML_biomarker_report();
	print gen_HTML_drug_sens_report();
	print gen_HTML_biomarker_evidence_report();
	print gen_HTML_therapy_summary_report();
	print gen_HTML_preferential_trial_report();
	print <<HTMLTAIL
</body></html>
HTMLTAIL
;

} elsif ($f_from_web) {
	print gen_HTML_biomarker_report();
	print gen_HTML_drug_sens_report();
	print gen_HTML_biomarker_evidence_report();
	print gen_HTML_therapy_summary_report();
	print gen_HTML_preferential_trial_report();
	print "<div class=debug>\n";
	print "<h2>Debug output:</h2>\n";
	print $Pottr->{'ruleset'}->{'debug_output'} ;
	print "</div>\n";

	print "<div class=debug>\n";
	print "<h2>POTTR Debug output:</h2>\n";
	print mk_trial_href( $Pottr->{'debug_output'} );
	print "</div>\n";
} elsif ($f_tsv) {
	print gen_biomarker_report() if $f_interp_variants ; 
	print gen_drug_sens_report() if $f_list_therapies;
	print gen_biomarker_evidence_report() if $f_list_evidence ; 
	print gen_therapy_summary_report() if $f_list_therapy_summaries ; 
	print gen_preferential_trial_report() if $f_list_trials; 
} elsif (! $f_verbose) {

# 	print hl
# 	print "\e[1;37mINPUT\e[0m\n";
# 	print map {"$_\n" } @input;
# 	print "\n\n";
	
	print hl(37, $POTTR_title)."\n";
	print "$POTTR_subtitle\n";
	print "$POTTR_disclaimer\n\n";

	if ( $f_clinical_report ) {
	}
	
	print gen_catypes_terminal()."\n" if $f_interp_catype ;
	print gen_prior_therapies_terminal()."\n" if $f_interp_catype ;
	print gen_biomarker_terminal()."\n" if $f_interp_variants ; 
	print gen_drug_sens_terminal()."\n" if $f_list_therapies;
	print gen_biomarker_evidence_terminal()."\n" if $f_list_evidence ; 
	print gen_therapy_summary_terminal()."\n" if $f_list_therapy_summaries;
	print gen_preferential_trial_terminal()."\n" if $f_list_trials; 
}





#######################################################################################################################################
exit 0;


__END__


=pod

=head1 NAME

pottr - Precision Oncology Therapy and Trial Recommender

=head1 SYNOPSIS


 pottr [options] [catype:``cancer-type''] [molecular_alterations]

For example:

 pottr [options] 'catype:Breast Cancer' 'ESR1:protein_expression' 'TP53:R175H'
 pottr [options] 'catype:Colorectal Cancer' 'BRAF:V600E' 'TP53:R175H'
 pottr [options] 'catype:Non-small cell lung cancer' 'EGFR:L858R' 'EGFR:T790M'
  
 pottr [options]  # Reads cancer type (catype:) and molecular alterations from STDIN


=head1 DESCRIPTION

Precision Oncology Treatment and Trial Recommender (POTTR) is an open-source, Perl-based clinical decision support system (CDSS) for recommending treatment for patients with advanced cancers, based on the results of tumour molecular profiling. It is intended to produce outputs that can be used by oncology physicians to aggregate available evidence about potential treatment options.

POTTR integrates patient's genomic testing results to produce evidence-based recommendations. Based on molecular alterations, this tool:
(1) Ranks potential therapies based on a selected evidence tiering (level of evidence) system;
(2) Predicts drug sensitivity/resistance, taking into account all available evidence and variant interactions, and;
(3) Preferentially selects and ranks therapies and selects clinical trials based on biomarker-matched evidence for a given tiering system.

As an expert system, POTTR requires a number of ontologies and knowledge bases to operate. 
Before use, POTTR needs to be properly be configured (see CONFIGURATION section).

Note, POTTR does not implement a bioinformatics module to determine pathogenicity of individual variants by default. 
The oncogenicity of each variant should be determined separately before calling POTTR. 
In addition, POTTR is yet capable of handling complex clinical variables, apart from limited histology type determination. 
Full interpretation of clinical data is planned in future versions.

Internally, POTTR is built on a purpose-built expert-system shell, using production rules to make inference on genomic alteration, clinical characteristics, and clinical evidence. 
The core reasoning mechanism behind POTTR is an inference tracking engine that forms the basis of evidence assessment.


=head1 OPTIONS

Basic Options:

   -?    --help                   This help message 
         --man                    Full documentation
   -c    --config-file            Specifiy the configuration file
  
Options for formatting output :                               

   -w    --web-output             Produce HTML outputs, including debugging reports
         --webpage                Produce a HTML web page
   -p    --clinical-report        Produce a clinical report based on biomarker/molecular profile
         --title                  Title of the report (e.g., patient name)
   -t    --export-tsv             Export recommendations as in Tab-Separated Value (TSV) format 
                               
Options for displaying individual sections:                               

   -C    --interpret-catypes      Generation interpretation of cancer types after ontology expansion
   -V    --interpret-variants     Generation interpretation of gene variants after ontology expansion
   -E    --list-evidence          Generation list of evidence matched to the patient's tumour genomic profile
   -R    --list-therapies         Generation list of therapies matching patient's tumour genomic profile
   -S    --list-therapy-summary   Generation list of therapies summary (recommendation) matching patient's 
                                  tumour genomic profile
   -T    --list-trials            Generation list of clinical trials patient's tumour genomic profile

Options for debugging:                              

   -q    --quiet                  Quiet output (default)
   -v    --verbose  --debug       List all facts and reasoning iterations for debugging
   -r    --print-rules            List all production rules 
         --modules                Select specific POTTR modules to run (internal use only)
   
                               


=head1 INVOCATION

Data input into POTTR is either be passed through command line argument or via standard input (one line per alteration or attributes, separated by LF).

The format is ``catype:cancer_type'' or ``biomarker:alteration''. 
The format of biomarker and molecular alteration are determined by evidence knowledge base. 

Example of inputs that are recognised by POTTR:

=over 2

=item * 
catype:Breast Cancer                - specifies a specific cancer type

=item * 
catype:Solid tumour                 - specifies a solid tumour

=item * 
EGFR:L858R                          - a simple mutation (p. notation by default)

=item * 
BRCA1:oncogenic_mutation            -- an unspecified oncogenic mutation

=item * 
BRCA1:oncogenic_mutation,germline   -- as above but germline

=item * 
FGFR1:amplification                 -- an gene amplification event 

=item * 
CDKN2A:deletion                     -- a (homozygous) gene deletion event 

=item * 
ROS1:fusion                         -- a unspecified fusion 

=item * 
ALK:ALK--EML4 fusion                -- a specific fusion    

=item * 
tumuor_mutational_burden:high       -- a specific genomic feature as biomarker

=item * 
ESR1:protein_expression             -- a protein expression (ER IHC-positive)

=item * 
ERBB2:overexpression                -- a HER2 overexpression (HER2 IHC 3+)

=item * 
RB1:loss_of_protein_expression      -- a protein expression (RB1 IHC-negative)

=back 

=head1 OUTPUT

=over 

=item * Section: ``CANCER TYPES MATCHED'' - reports details of biomarker after ontology expansion.

=item * Section: ``BIOMARKER ALTERATIONS'' - reports details of biomarker after ontology expansion.

=item * Section: ``THERAPY SENSITIVITY/RESISTANCE'' - given a genomic profile and cancer type, summarises the report of sensitivity and resistance to therapy/therapy classes.

=item * Section: ``EVIDENCE OF BIOMARKER-MATCHED THERAPIES'' - list of evidences matched in knowledge base with respect to given biomarker alterations.

=item * Section: ``BIOMARKER-MATCHED CLINICAL TRIALS''

=over 2

=item * Rank 

=item * Trial ID - unique identifier of the trial 

=item * Drugs - Full title of the trial                                                                                

=item * Drugs - list of drugs (and drug combinations) studied in this trial                                  

=item * Drug classes  - list of drug classes (and combinations) studied in this trial

=item * Transitive class efficacy tier, Transitive efficacy tier, Drug maturity tier, Drug class maturity tier, Combo maturity tier, and Combo class maturity tier - metrics used for ranking preferentially identified trials.

=item * Health conditions  - matched cancer types

=item * Location / Postcodes                        

=back

=back

=head1 CONFIGURATION

The default configuration file is pottr_config. The path of configuration file can be change 
by using the ``-c'' option via the command line.

An example configuration file is listed in pottr_config.

=head2 GLOBAL OPTIONS [REQUIRED] 

Synopsis:
  
  data_dir  =  data/
  cache_dir =  data/

Specification of ``data'' directories (database files) and ``cache'' directories (for storing intermediate files). 
The databases files are relative to these directories specified here.

=head2 DISEASE ONTOLOGY [REQUIRED] 

Synopsis:
  
  disease-ontology-file = doid.obo

Disease ontology file is in the OBO ontology format, and the latest version can be obtained from https://disease-ontology.org/. 
Please follow the licesing requirements for details.
 
Note future versions will use a dedicated ontology for precision oncology to more closely align with clinical practice.

=head2 VARIANT FEATURE FILE [OPTIONAL] 

Synopsis:
  
  variant-feature-file = feature_lookup.txt

This database is used for annotating gene features clinically relevant to decision making.

File format: tab-separated values (TSV) without header comprising entries:

=over 2

=item  * gene     - (string) HGNC gene symbol

=item  * feature  - (string) gene feature, e.g., exon or specific protein domains.

=item  * start_aa - (integer) starting aa number.

=item  * stop_aa  - (integer) stopping aa number.

=item  * comments - (string) user-specified comments. not processed 

=back

  Example entry:
  ABL1    exon_7  382     443
    

=head2 ONCOGENICITY INFERENCE/ASSERTION FILE [OPTIONAL]
 
Synopsis: 
  
  oncogenicity-rules-file = variant_annotation_rules.txt 

This optional database contains rules used for optional variant interpretation. This 
step is not essential for decision making if pathogenicity of a variant has been determined 
before POTTR is called. 

The following example demonstrates an excerpt translated from OncoKB pathogenicity database:

  EGFR:L858R => EGFR:oncogenic_mutation [CERTAIN:oncogenicity]
  EGFR:L858R => EGFR:gain-of-function_mutation [CERTAIN:mutation_effect]


=head2 BIOMARKER EVIDENCE FILE [REQIURED] 

Synopsis: 
  
  therapy-evidence-database-name = TOPOGRAPH
  therapy-evidence-database-file = GIMRKB-master.tsv

This section specifies core evidence database. Multiple sections of evidence files are allowed.
However, the tiering system must be identical in all databases. 
Evidence from all knowledge bases will be assimilated and assessed.

FORMAT: tab-separated values (TSV) 

Fields:

=over 2

=item * Date        - (String) Date of the particular entry, used for versioning.

=item * Tier        - (Requied) The tier of a particular evidence entry (See below)

=item * Biomarker   - Name of the biomarker, such genes (e.g., EGFR), protein (e.g., 
              CD274) or special biomarkers (e.g., Microsatellite instability).
              Multiple alterations are allowed and separated by "+" (WITHOUT 
              whitespace)

=item * Alteration  - (String) gene mutation (using HGVS nomenclature, p. notation by 
              default if not specified), alteration of protein expression (e.g., 
              ERBB2 overexpression), fusion (e.g., BCR-ABL1 fusion), etc. Multi-
              ple alterations are allowed and separated by comma or semicolon. 
              Co-alteration (of the same biomarker/gene) are joint by " and ".
              Co-alteration (of the different biomarker/gene) are also joint 
              by " and ", but will preceed with. Negation operators such as 
              "except" or "not" are also allowed (e.g., "BRAF:alteration and 
              NOT BRAF:V600")

=item * Tumour Type - (Also allowed: Cancer type). Specification of cancer type 
              matching evidence Multiple tumour types allowed, separated by ";"

=item *  Drugs       - List of therapies separated buy ";". Drug combinations should be
               separated by " + "

=item *  Comments    - (string) Not processed

=item *  Evidence    - (string) List of evidence separated by ", ". Suggested using PMID 
              and DOI (for publications) or NCT number (for early phase trials)

=back


=head2 EVIDENCE GRADING SECTION [REQUIRED]

Defining tiers for ranking treatment and clinical trial. The following section defines an example of configuration using OncoKB-like tiering system (i.e., 1 2 3 4 R1 R2).

  set tiers-list:1 1R 1B 1BR 2 2R 3 3R 3B 3BR 4 4R 4B 4BR R1 R2 R2B # Specification/list of tiers
  set tiers-list-resistance:R1 R2 R2B  # Specification/list of drug resistance tiers

Custom specification of different tier ranking methods: 

  set tier-ranking-mode:pedantic        => tier-rank-order:R1 R2 R2B 1 1R 1B 1BR 2 2R 3 3R 3B 3BR 4 4R 4B 4BR
  set tier-ranking-mode:conservative    => tier-rank-order:R1 1 R2 R2B 1R 1B 1BR 2 2R 3 3R 3B 3BR 4 4R 4B 4BR
  set tier-ranking-mode:liberal         => tier-rank-order:R1 1 1R 1B 1BR 2 2R 3 R2 R2B 3R 3B 3BR 4 4R 4B 4BR


=head2 THERAPY DATABASE [REQUIRED] 

Synopsis:
  
  therapy-database-file = drug_database.txt                   
  
The therapy database consists of a TSV file with two fields:

=over 2 

=item * drug        - (string) preferred drug name, followed by synonyms, separated by '|'

=item * drug_class  - (string) a string showing the immediate (primary) drug class name. 
  Multiple drug classes are allowed, separated by '; '

=back
    
Example: 
  
  Afatinib|BIBW 2992|BIBW2992	EGFR_inhibitor,second_generation
  Alectinib|Alecensa|AF-802|AF802|CH5424802|RG7853|RO5424802	ALK_inhibitor,second_generation
  Savolitinib|AZD 6094|AZD6094|HMPL-504|Volitinib	MET_inhibitor,type_1
  Entrectinib|RXDX 101|RXDX-101|RXDX101|Rozlytrek ALK_inhibitor,third_generation; TRK_inhibitor,first_generation; ROS1_inhibitor


=head2 THERAPY HIERARCHY DATABASE [REQUIRED] 

Synopsis:
  
  drug-class-hierarchy-file = drug_class_hierarchy.txt            


A simple, line-based, tab-separated database of taxonomic hierarchy of of therapy classes (is-a relationship). 
This is required for both therapy and trial search/matching.
 
  For example, the following entry:
  immune_checkpoint_blockade  <TAB>  immune_checkpoint_blockade,PD-1_targeting  <TAB>   anti-PD-1_monoclonal_antibody

indicate

  anti-PD-1_monoclonal_antibody IS-A immune_checkpoint_blockade,PD-1_targeting
  
and 

  immune_checkpoint_blockade,PD-1_targeting IS-A immune_checkpoint_blockade      

respectively



=head2 CLINCAL TRIALS DATABASE [OPTIONAL]

Synopsis:
  
  clinical-trial-database-file = ANZCTR-potential-trial-summary.tsv  

Specify this file (or files) if clinical trial matching is required. 

The database (in TSV format) with title consists of the following fields:

=over 2

=item * trial_id                     -  (string) Unique Clinical trial ID

=item * trialacronym                 -  (Optional) clinical trial acronym

=item * phase                        -  (string) Phase of clinical trial

=item * studytitle                   -  (string) Title of the study

=item * drug_list                    -  (strings) List of study drugs for matching, 
                                   separated by semicolons. Drug combinations 
                                   are allowed, separated by ' + ' (NB: spaces 
                                   before and after '+' sign are required)

=item * drug_classes                 -  (strings) List of drug classes (of the study 
                                   drugs), separated by semicolon. Drug classes 
                                   are defined in therapy database. Drug class 
                                   combinations are allowed.
=item * combo_list                   -  (strings) List of study drug combinations 
                                   specified in the trial for matching, separated 
                                   by semicolons. Drug combinations are separated by 
                                   ' + ' (NB: spaces before and after '+' sign are required)
                                   
=item * combo_classes                -  (strings) List of drug classes (of the study 
                                   drug combinations), separated by semicolon. The specific 
                                   drug classes must be defined in the therapy database. 
=item * healthcondition (aka catype) -  (strings) List of cancer types for matching, 
                                   separated by semicolons. 

=item * recruitmentstatus -  (string) Recruitment status

=item * postcode -  (string) Location/ZIP/Postcode. 

=item * ext_weblink -  (string:URL) Link to external website (e.g., 
                                    clinicaltrials.gov or ANZCTR web page)

=back


=head2 TRIAL ELIGIBITIY FILTERING FILE [OPTIONAL]

Synopsis: 

  clinical-trial-eligibility-file = ANZCTR-trials-eligibility.tsv


TSV file specifying criteria specifying:

=over 2

=item * trial_id               - (string) Unique Clinical trial ID

=item * eligibility_criteria   - trial filtering using POTTR rules

=back 

  For example:
  trial_id	eligibility_criteria
  NCT02857270	NRAS:alteration; catype:Melanoma
  NCT02857270	catype:Non-small cell lung cancer; (BRAF:alteration OR KRAS:alteration OR NRAS:alteration)
  NCT02857270	catype:colorectal cancer; BRAF:V600E


=head2 THERAPY MATURITY FILE - [OPTIONAL] (EXPERIMENTAL) 

  Synopsis:
  
  therapy-maturity-database-file = drug_trial_phases_CTGOV.txt

Databases used to assess the maturity of therapy/therapy classes with respect to the stage of development.
These databases use TSV format (without title) containing the following fields in order:

=over 2

=item * trial_id  - unique trial identified or e.g., NCT0XXXXXXX or ACTRNXXXXXXXXX

=item * phase     - Phase of the trial (e.g., Phase 1, Phase 2, Phase 3, Phase 4)

=item * attribute - one of the following
               "drug" - drug listed in the trial
               "drug_class" - drug classes of the drug listed in the trial
               "combo" - combination of drugs used in the trial
               "combo_class" - classes of combination of drugs used in the trial

=item * therapy   - for drug (or drug class) combinations, separate each ingredient by 
               either " + " or "|"

=back 

  Example:
  NCT02151981     Phase 3 drug    Osimertinib
  NCT02151981     Phase 3 drug    Pemetrexed
  NCT02151981     Phase 3 drug_class      EGFR_inhibitor,third_generation
  NCT02151981     Phase 3 drug_class      antimetabolite
  NCT02151981     Phase 3 drug_class      platinum-based_antineoplastic_agent
  NCT02151981     Phase 3 combo   Carboplatin|Cisplatin|Pemetrexed
  NCT02151981     Phase 3 combo   Osimertinib
  NCT02151981     Phase 3 combo_class     antimetabolite|platinum-based_antineoplastic_agent
  NCT02151981     Phase 3 combo_class     EGFR_inhibitor,third_generation


=head1 RETURN VALUE

Returns 0 on success, non-zero if error


=head1 DIAGNOSTICS

For debugging, use ``-v'' switch to examine the rules activated and facts at each iteration.


=head1 EXAMPLES

=over

=back

=head1 CITING POTTR


=head1 DISCLAIMER AND WARNINGS

Specialist knowledge in oncology is constantly evolving. Making treatment recommendations  
for cancer patients requires  careful evaluation  of clinical situations in its entirety.  
Results produced by POTTR are STRICTLY FOR RESEARCH USE ONLY  and  CANNNOT SUBSITUTE  for 
professional oncology advice. Using POTTR is WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND. 
Inappropriate use of this software may lead to harm.


=head1 AUTHOR

Frank Lin. MBChB PhD FRACP. Kinghorn Centre for Clinical Genomics, Garvan Institute of Medical Research, Sydney, Australia.
For bug fixes, please contact: f dot lin at garvan dot org dot au. 

=head1 COPYRIGHT AND LICENSE

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 SEE ALSO

=cut
