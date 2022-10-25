#!/usr/bin/perl
##############################################################################
#
# CancerTypes.pm - precision oncology trial and therapy recommender
# 
# Cancer type ontology - a simple cancer type ontology
#
# This file is part of Oncology Treatment and Trials Recommender (OTTR) and 
# Precision OTTR (POTTR).
#
# Copyright 2022, Frank Lin & Kinghorn Centre for Clinical Genomics, 
#                 Garvan Institute of Medical Research, Sydney
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

package CancerTypes;

use strict;
use warnings;

use POSIX;

use lib '.';
use TSV;

use re 'eval';
use POTTRConfig;
use Storable;

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


sub trim { my $s = shift; $s =~ s/^\s*|\s*$//g if (defined $s) ; return $s; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }

our $f_initiailised = undef;

our $catype_regex_M ;
our %catype_preferred_term;
our %catype_synonyms;
our %catype_signature_cache;
our %catype_is_a;
our %catype_has_child;
our %catype_not_found;
our %preferred_catype_term_cache;

sub mk_signature {
	my $x = shift;
	return $catype_signature_cache{$x} if exists $catype_signature_cache{$x} ;
	$x =~ s/[^A-Za-z0-9]//g;
	$x = lc($x);
	$catype_signature_cache{$x} = lc($x);
	return lc($x);
}

sub get_preferred_catype_term {
	&ON_DEMAND_INIT;
	
	my $x = shift;
	
	return $preferred_catype_term_cache{$x} if exists $preferred_catype_term_cache{$x} ;
	
	if ( exists $catype_preferred_term{ (my $xsig = mk_signature($x)) } ) {
		return ( $preferred_catype_term_cache{$x} = $catype_preferred_term{ $xsig } ) 
	}
	
	if ( defined $x ) {
		warn "\e[1;31m".(my $msg = "CancerTypes.pm: ``$x'' not found in the database.")."\e[0m\n";
		main::error_hook( $msg ) if ( UNIVERSAL::can("main", "error_hook") ); # Generating an alert
		$catype_not_found{ $x } ++ ;
		$catype_preferred_term{ mk_signature($x) } = $preferred_catype_term_cache{$x} = $x;
	}
	return $x;
}


sub gen_stats {
	my %seen;
	my %parents;
	my %children_list;
	my %offsprings_list;
	my %siblings_list;
	
	for my $a (sort keys %catype_is_a) {
		my @b = keys %{ $catype_is_a{$a} };
		$seen{$_} = 1 for @b;
		$seen{$a} = 1;
		$parents{$a} = scalar( @b );
		for my $b (@b) {
			$children_list{$b}{$a} = 1;
		}
	}
	for my $a (sort keys %seen) {
		my @ancestors = get_ancestors( get_preferred_catype_term($a) );
		$offsprings_list{ mk_signature($_) }{$a} = 1 for @ancestors ;
		if ( exists $catype_is_a{$a} ) {
			my @p = keys %{ $catype_is_a{$a} };
			for my $p (@p) {
				for my $s ( keys %{ $children_list{$p} } ) {
					$siblings_list{$a}{$s} = 1;
				}
			}
		}
	}
	
	my %children = map { $_ => scalar keys %{ $children_list{$_} } } keys %children_list;
	my %offsprings = map { $_ => scalar keys %{ $offsprings_list {$_} } } keys %offsprings_list ;
	my %siblings = map { $_ => scalar keys %{ $siblings_list {$_} } } keys %siblings_list ;
	my %children_siblings_ratio = map { $_ => ( ($children{$_} // 0) / ($siblings{$_} // 1) )  } keys %seen;
	my %offsprings_siblings_ratio = map { $_ => ( ($offsprings{$_} // 0) / ($siblings{$_} // 1) )  } keys %seen;
	
# 	for my $a (sort { ($children_siblings_ratio{$b} // 0) <=> ($children_siblings_ratio{$a} // 0) || $a cmp $b } keys %seen) {
	for my $a (sort { ($offsprings_siblings_ratio{$b} // 0) <=> ($offsprings_siblings_ratio{$a} // 0) || $a cmp $b } keys %seen) {
# 		print STDERR join("\t", $parents{$a} // 0, $siblings{$a} // 0, $children{$a} // 0, $offsprings{$a} // 0, sprintf( "%.3f", $offsprings_siblings_ratio{$a} // 0), get_preferred_catype_term($a) )."\n";
	}
}

sub load_catype_database($) {
	my $srcf = shift;
	for my $line ( file($srcf) ) {
		chomp $line ;
		my @parts = split /\t/, $line ;
		my $prev_catype = undef;
		for my $p (@parts) {
			my ($term, @synonym) = map { trim($_) } split /\|/, $p;
			next if $term =~ /^\s*$/;
			for my $t ($term, @synonym) {
				next if $t =~ /^\s*$/;
				$catype_synonyms{ mk_signature($term) }{ $t } = 1;
# 				print "$term (".mk_signature($term).") = $t\n";
				$catype_preferred_term{ $t } = $term;
				$catype_preferred_term{ mk_signature($t) } = $term;
			}
			if ( defined $prev_catype ) {
				$catype_is_a{ mk_signature($term) }{ mk_signature($prev_catype) } = 1 ;
				$catype_has_child{ mk_signature($prev_catype) }{ mk_signature($term) } = 1 ;
			}
			$prev_catype = $term;
		}
	}
	
# 	gen_stats;
}

sub is_a { # catype is-a catype
	&ON_DEMAND_INIT;
	my $x = mk_signature($_[0]);
	my $y = mk_signature($_[1]);
	my $visited = $_[2];
	my %visited_stack;
	
	return 1 if $x eq $y;
	
	$visited = \%visited_stack if ! defined $visited ;
	${$visited}{$x} = 1;
	
# 	print "$x\t$y\n";
	if ( exists $catype_is_a{$x} ) {
		return 1 if exists $catype_is_a{$x}{$y};
		for my $z ( keys %{ $catype_is_a{$x} } ) {
			next if exists ${$visited}{$z};
			return 1 if is_a($z, $y, $visited);
		}
	}
	return 0;
}


sub get_ancestors { 
	&ON_DEMAND_INIT;
	
	my $x = mk_signature( get_preferred_catype_term($_[0]) );
	my $visited = $_[1];
	my @visited_stack;
	$visited = \@visited_stack if ! defined $visited ;
	push @{ $visited }, $x;
	
# 	print "\e[1;31m$x\e[0m\n";
	if ( exists $catype_is_a{$x} ) {
# 		print STDERR "!";
# 		print "\e[1;31m$x\t$catype_is_a{$x}\e[0m\n";
		for my $z ( keys %{ $catype_is_a{$x} } ) {
			next if ! defined $z;
			next if grep { $z eq $_ } @{$visited};
# 			print "\e[1;31m\t$z\e[0m\n";
			get_ancestors($z, $visited );
		}
	}
# 	print "\e[1;31m".join(" | ", ( map { "$_ (". $catype_preferred_term{$_}. ")" } keys %visited_stack ) )."\e[0m\n";
	my @a = map { $catype_preferred_term{$_} } @visited_stack ;
	return @a;
}



sub get_offsprings { 
	&ON_DEMAND_INIT;
	
	my $x = mk_signature( get_preferred_catype_term($_[0]) );
	my $visited = $_[1];
	my @visited_stack;
	$visited = \@visited_stack if ! defined $visited ;
	push @{ $visited }, $x;
	
	if ( exists $catype_has_child{$x} ) {
		for my $z ( keys %{ $catype_has_child{$x} } ) {
			next if ! defined $z;
			next if grep { $z eq $_ } @{$visited};
			get_offsprings($z, $visited );
		}
	}
	my @a = map { $catype_preferred_term{$_} } @visited_stack ;
	return @a;
}


sub get_all_catypes { 
	&ON_DEMAND_INIT;
	return sort map { $catype_preferred_term{$_} } keys %catype_preferred_term;
}

our $sep_regex = "[[:space:]\\-]?";

sub get_catype_regex {
	&ON_DEMAND_INIT;
	my @elems ;
	
	for my $arg (@_) {
		push @elems, keys %{ $catype_synonyms{ mk_signature( get_preferred_catype_term($arg) ) } } ;
	}
	
	my %by_signature;
	for my $e (@elems) {
		$by_signature{mk_signature($e)}{$e} ++;
	}
	
	@elems = ();
	for my $k (sort keys %by_signature) {
		my @a = keys %{ $by_signature{$k} };
		my @m = @a;
		if (scalar @a >= 2 ) {
			my ($s) = sort { length $b <=> length $a } grep { /[\- ]/ } @a;
# 			warn $s;
			if ( defined $s ) {
				$s =~ s/[ \-]/$sep_regex/g;
				my $d = grep { /^$s$/i } @a ;
				if ( $d == scalar(@a) ) {
					@m = ($s);
				}
			}
		}
		push @elems, @m;
	}

	my $regex = join("|", sort { length $b <=> length $a } @elems );

	return $regex ;
}


sub catype_regex_escape { return $_[0] =~ s/(['])/\\$1/r; }

my @word_synonym_set = (
	'(?:lung|pulmonary)',
	'(?:stomach|gastric)',
	'(?:kidney|renal)',
	'(?:liver cell|hepatocellular)',
	'(?:uterine|uterus|endometrial|endometrium)',
	'(?:harbouring|with)',
	' of(?: the)? ',
);

sub get_catype_regex_M_expand_regex {
	my $x = $_[0];
	$x =~ s/(?:cancer|tumou?r)/(?i:malignant )?(?i:cancer|tumou?r|carcinoma|sarcoma|neoplasm|malignancy|malignancies)s?/i;
	for my $r (@word_synonym_set) {
		$x =~ s/$r/$r/i;
	}
	$x =~ s/(appendi|pancrea|prostat|uter|cervi|endometri|ovar|oesophag|lymph|gynaecologic|diffuse|metastas|mutat)(?:s|e|[ie]s|t?ic|ine|ion|oid|ed|ant|oma|us|c?e?al|x|um|ian|y|d)/$1(?:s|e|[ie]s|t?ic|ine|ion|oid|ed|ant|oma|us|c?e?al|x|um|ian|y|d)/i;
	$x =~ s/(o)(eso)/$1?$2/i;
	$x =~ s/(o)(u)/$1$2?/i;
	$x =~ s/([ao])(e)/$1?$2/i;
	$x =~ s/(c)(y)$/$1(?:$2|ie)/i;
	$x =~ s/\b(?:p\.)/(?:p\.)?/;
	$x =~ s/[, \-]+/.{0,2}/g;
	$x =~ s/,\s*NOS$//;
	$x .= 's?' if $x !~ /s\??$/i;
	return $x;
}

sub get_catype_regex_M {
	&ON_DEMAND_INIT;
	my @elems ;
	for my $arg (@_) {
		push @elems, keys %{ $catype_synonyms{ mk_signature( get_preferred_catype_term($arg) ) } };
	}
	my %by_signature;
	for my $e (@elems) {
		$by_signature{mk_signature($e)}{$e} ++;
	}
	
	@elems = ();
	my %pref_name ;
	for my $k (sort keys %by_signature) {
		my @a = keys %{ $by_signature{$k} };
		my @m = @a;
		$pref_name{$_} = get_preferred_catype_term($_) for @m ;
		if (scalar @a >= 2 ) {
			my ($s) = sort { length $b <=> length $a } grep { /[\- ]/ } @a;
			if ( defined ($s) ) {
				$s =~ s/[ \-]/.?/g; #$sep_regex
				my $d = grep { /^$s$/i } @a ;
				if ( $d == scalar(@a) ) {
					@m = ($s);
					$pref_name{$s} = get_preferred_catype_term($a[0]);
				} 
			}
		} 
		push @elems, @m;
	}

	my $regex = join("|", map { get_catype_regex_M_expand_regex($_)."(?{\$MMM='".catype_regex_escape($pref_name{$_})."'})" } sort { length $b <=> length $a } grep { length } @elems);
	
# 	print STDERR "".($regex  =~ s/\)\|/\)\n\t|/gr)."\n";
	
	return $regex ;
}


sub get_catype_synonyms {
	&ON_DEMAND_INIT;
	return sort keys %{ $catype_synonyms{ mk_signature( get_preferred_catype_term($_[0]) ) } };
}

our %match_catype_whole_word_cache;
our %match_catype_cache;

our $match_catype_whole_word_cache_fn ; # = POTTRConfig::mk_type_path('cache', 'pottr_catypes_whole_word.cache');
our $match_catype_cache_fn ; # = POTTRConfig::mk_type_path('cache', 'pottr_catypes.cache');

sub match_catype_whole_word {
	&ON_DEMAND_INIT;
	my $string = $_[0];
	my $strict = $_[1]; # if strict flag is set, return undef.
	
# 	print STDERR "CancerTypes::match_catype_whole_word: $string\n";
# 	$catype_regex = get_catype_regex( get_all_catypes() )  if ! defined $catype_regex ;
# 	return get_preferred_catype_term($x) ; # if $string =~ /^\s*($catype_regex)\s*$/i ;
	return $match_catype_whole_word_cache{$string} if exists $match_catype_whole_word_cache{$string};
	
	$catype_regex_M = get_catype_regex_M( get_all_catypes() )  if ! defined $catype_regex_M ;

	if ( -f $match_catype_whole_word_cache_fn ) {
		%match_catype_whole_word_cache = map { chomp; my ($a, $b) = split /\t/, $_; $a => $b } file($match_catype_whole_word_cache_fn); 
	}

	our $MMM = '';
	if ( $string =~ /^[[:punct:][:space:]]*(${catype_regex_M})[[:punct:][:space:]]*\s*$/ig ) {
		my $m = $match_catype_whole_word_cache{$string} = get_preferred_catype_term($MMM);
# 		print STDERR "\tCancerTypes::match_catype_whole_word: $string: new: $1: $MMM => $m\n" ;
		open FOUT_CACHE, ">$match_catype_whole_word_cache_fn";
		print FOUT_CACHE map { "$_\t$match_catype_whole_word_cache{$_}\n" } grep { defined $match_catype_whole_word_cache{$_} } (sort keys %match_catype_whole_word_cache ); # die $match_catype_whole_word_cache{$_} if ! exists $match_catype_whole_word_cache{$_} 
		close FOUT_CACHE;
		return ($m);
	}
	
	return undef if defined $strict;

	return ($match_catype_whole_word_cache{$string} = get_preferred_catype_term($string) );
}



sub match_catype {
	&ON_DEMAND_INIT;
	my $string = $_[0];
	my %x ;
	
	return %{ $match_catype_cache{$string} } if exists $match_catype_cache{$string};

	if ( -f $match_catype_cache_fn ) {
		%match_catype_cache = %{ retrieve $match_catype_cache_fn };
	}
	
	$catype_regex_M = get_catype_regex_M( get_all_catypes() )  if ! defined $catype_regex_M ;
# 	$catype_regex_M = get_catype_regex_M_expand_regex( get_all_catypes() )  if ! defined $catype_regex_M ;
	
	our $MMM = '';
	while ( $string =~ /(${catype_regex_M})/ig ) {
		my $start = length($`);
		my %match = (start => $start, end => $start + length($1), len =>length($1), text => $1 );
		my $term = get_preferred_catype_term($MMM);
		push @{ $x{$term} }, \%match;
	}
	
	my %wc = %x;
	
	$match_catype_cache{$string} = \%wc;

	store \%match_catype_cache, $match_catype_cache_fn;
	
	return %x;
}


sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;
	
	$match_catype_whole_word_cache_fn = POTTRConfig::mk_type_path('cache', 'pottr_catypes_whole_word.cache');
	$match_catype_cache_fn = POTTRConfig::mk_type_path('cache', 'pottr_catypes.cache');
	
	for my $db ( POTTRConfig::get_paths('data', 'cancer-types-ontology-file') ) {
		if ( ! -f $db ) {
			warn "\e[1;31mCancerTypes.pm: ONTOLOGY WARNING: $db not found.\e[0m\n" ;
			next;
		}
		load_catype_database $db;
	}
}

# 
# load_cancer_type_databases('data/cancer_types.txt');
# 
# print get_preferred_catype_term('Malignant melanoma')."\n";
# print join("; ", get_ancestors('Well-differentiated liposarcoma'))."\n";
# print is_a('Well-differentiated liposarcoma', 'Solid Tumour')."\n";
# print is_a('Well-differentiated liposarcoma', 'Liposarcoma')."\n";
# print is_a('Well-differentiated liposarcoma', 'Breast cancer')."\n";
# print is_a('Well-differentiated liposarcoma', 'Well-differentiated liposarcoma')."\n";
# 
# print get_all_catypes();
# print "\n";
# # print ($regex_M =~ s/\|/\n\t\|/gr)."\n";
# 
# 
# 
# my $string  = "Aggressive fibromatosis is also known as Desmoid, amongst rare types of sarcomas. Triple-negative breast cancer is a type of breast cancer. Non-Hodgkin's lymphoma is a haematologic malignancy, not a solid tumor.";
# 
# my %results = match_catype($string);
# 
# print "$string\n";
# for my $catype_term (sort keys %results) {
# 	for my $match ( @{ $results{$catype_term } } ) {
# 		print "\e[1;36m$catype_term\e[0m => \e[1;35m($$match{start}, $$match{end}, '$$match{text}')\e[0m\n";
# 	}
# }
# 
# for my $test_word (
# 	"Liposarcomas",
# 	"Solid tumours",
# 	"Acute myeloid leukemias",
# 	"Glioblastoma multiforme",
# 	"Basal cell carcinoma",
# 	) {
# 	my $m = ( match_catype_whole_word($test_word) // "Unmatched") ;
# 	print join(" => ", "\e[1;37m$m\e[0m", get_ancestors($m) )."\n";
# }

# 
# my $regex_M = get_catype_regex_M( get_all_catypes() );
# my $M;
# while ( $string =~ /($regex_M)/ig ) {
# 	print "\e[1;36m$1 -> ".get_preferred_catype_term($1)."\e[0m: $M";
# 	print "\t\e[1;34m".join("; ", get_catype_synonyms($1))."\e[0m\n";
# }

1;
