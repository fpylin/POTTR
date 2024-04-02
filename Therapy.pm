#!/usr/bin/perl
##############################################################################
#
# Therapy.pm - precision oncology trial and therapy recommender
# 
# Therapy assessment functions
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

package Therapy;

use strict;
use warnings;

use POSIX;

use lib '.';
use POTTRConfig;

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

## Utility functions #################################
sub min { return undef if (! scalar(@_) ); my $v = shift; for (@_) { $v = $_ if ($_ < $v) ; } return $v; }
sub max { return undef if (! scalar(@_) ); my $v = shift; for (@_) { $v = $_ if ($_ > $v) ; } return $v; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }
sub trim { my $s = shift; $s =~ s/^\s*//; $s =~ s/\s*$//; return $s; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub levenshtein {
	my ($s1, $s2) = @_;
	my ($len1, $len2) = (length $s1, length $s2);
	return $len2 if ($len1 == 0);
	return $len1 if ($len2 == 0);
	my %mat;
	for (my $i = 0; $i <= $len1; ++$i) { 
		for (my $j = 0; $j <= $len2; ++$j) {
			$mat{$i}{$j} = 0;
			$mat{0}{$j} = $j;
		}
		$mat{$i}{0} = $i;
	}
	my @ar1 = split(//, $s1);
	my @ar2 = split(//, $s2);
	for (my $i = 1; $i <= $len1; ++$i) {
		for (my $j = 1; $j <= $len2; ++$j) {
			my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
			$mat{$i}{$j} = min( $mat{$i-1}{$j} + 1, $mat{$i}{$j-1} + 1, $mat{$i-1}{$j-1} + $cost );
		}
	}
    return $mat{$len1}{$len2};
}

######################################################################
our $f_initiailised = undef;

our %known_drug_combinations;    # a register of known drug combination

our %drugs_not_found;

our %drug_class = () ;           # drug_signature -> drug_class_signature 
our %drug_class_name = () ;      # drug_class_signature -> drug_class_name
our %drug_class_is_a = () ;      # ''drug_signature''  is-a  ''drug_signature''
our %drug_class_has_parent = (); # ''drug_signature''  is-a  ''drug_signature''
our %drug_name = ();             # drug_name -> primary/preferred_drug_name
our %drug_preferred_name = () ;  # drug_signature -> drug_name
our %drug_synonyms = () ;        # drug_signature -> drug_name
our %drug_combo_regimens = () ;

our $drug_regex = undef;
our %drug_signature_cache = ();

our @hierarchy_pairs;

sub tidy_tier { my $tier = shift; return '4U' if $tier =~ /^[US]$/; return $tier ; }

######################################################################
sub mk_signature {
	my $x = shift;
	return $drug_signature_cache{$x} if exists $drug_signature_cache{$x} ;
	$x =~ s/[^A-Za-z0-9]//g;
	my $lc_x = lc($x);
	$drug_signature_cache{$x} = $lc_x;
	return $lc_x;
}

sub is_a_drug { 
	&ON_DEMAND_INIT;
	return exists $drug_name{ mk_signature($_[0]) } ; 
}

sub is_a_drug_class { 
	&ON_DEMAND_INIT;
	return exists $drug_class_name{ mk_signature($_[0]) } ; 
}

sub is_a { # drug is-a drug_class  
	&ON_DEMAND_INIT;
	my $x = mk_signature($_[0]);
	my $y = mk_signature($_[1]);
	my $visited = $_[2];
	my %visited_stack;
	
	$visited = \%visited_stack if ! defined $visited ;
	${$visited}{$x} = 1;
	
# 	print "$x\t$y\n";
	if ( exists $drug_class_has_parent{$x} ) {
		return 1 if exists $drug_class_has_parent{$x}{$y};
		for my $z ( keys %{ $drug_class_has_parent{$x} } ) {
			next if exists ${$visited}{$z};
			return 1 if is_a($z, $y, $visited );
		}
	} elsif ( exists $drug_class{$x} ) {
		return 1 if exists $drug_class{$x}{$y};
		for my $z ( keys %{ $drug_class{$x} } ) {
			next if exists ${$visited}{$z};
			return 1 if is_a($z, $y, $visited);
		}
	}
	return 0;
}

sub is_drug_PBS_reimbursed { 
	&ON_DEMAND_INIT;
	return is_a($_[0], 'PBS_listed') ; 
}

sub is_drug_eviQ_listed { 
	&ON_DEMAND_INIT;
	return is_a($_[0], 'eviQ_listed') ; 
} 


################################################################################
{package CombinationCounter;
        
sub new {
	my ($class, $groups_ref) = @_;
	my @a = map { 0 } ( 0 .. $#{ $groups_ref } );
	my $self = { 'groups' => $groups_ref, 'value' => \@a };
	bless($self, $class);
	return $self;
	}

sub incr {
	my $self = shift;
	my @a = @{ $self->{'value'} } ;
	my @g = @{ $self->{'groups'} } ;
	$a[0] ++;
	for (my $d=0; $d <=$#a; ++$d) {
		if ( ($d < scalar(@g)) and ($a[$d] >= scalar @{ $g[$d] }) ) {
			$a[$d] = 0;
			$a[$d+1] = 0 if $d >= scalar(@a);
			$a[$d+1] ++;
		}
	}
	$self->{'value'} = \@a;
}

sub overflow {
	my $self = shift;
	my @a = @{ $self->{'value'} } ;
	my @g = @{ $self->{'groups'} } ;
	return ( (scalar(@g)<scalar(@a)) and ($a[scalar(@g)] > 0) )
}

sub get {
	my $self = shift;
	my @a = @{ $self->{'value'} } ;
	my @g = @{ $self->{'groups'} } ;
	my @v = map { $g[$_][ $a[$_] ] } (0..$#g);
	return @v;
}

1;
}

sub get_all_parents { # drug
	&ON_DEMAND_INIT;
	
	my $x = mk_signature($_[0]);
	my $visited = $_[1];
	my %visited_stack;
# 	print STDERR "\e[1;31m$_[0]\e[0m\n";
	$visited = \%visited_stack if ! defined $visited ;
	${$visited}{$x} = 1;
	
	if ( exists $drug_class_has_parent{$x} ) {
		for my $z ( keys %{ $drug_class_has_parent{$x} } ) {
			next if exists ${$visited}{$z};
			get_all_parents($z, $visited );
		}
	} elsif ( exists $drug_class{$x} ) {
		for my $z ( keys %{ $drug_class{$x} } ) {
			next if exists ${$visited}{$z};
			get_all_parents($z, $visited);
		}
	}
	
# 	print ">> $_[0]\t$drug_class_name{$x}\n" if exists $drug_class_name{$x};

	if ( exists $drug_class_name{$x} and $drug_class_name{$x}  =~ / \+ / ) { # if the treatment is a combination therapy, then also traverse individual drug's parents
		my @components =  split / +\+ +/m, $drug_class_name{$x};
		my @components_r;
		for my $z ( @components ) {
			next if exists $$visited{$z};
			get_all_parents( $z, $visited);
			push @components_r, [ get_all_parents($z) ];
# 			print map {">>>> $_ @components_r\n"} get_drug_classes($z); 
		}
		
		for ( my $combo = CombinationCounter->new(\@components_r); ! $combo->overflow(); $combo->incr() ) {
			my @combo = $combo->get();
# 			print map { "[[$_]]\n" } @combo;
			my $z = join(" + ", ( sort map { get_normalised_treatment_class_name($_) } @combo) );
			my $zsig = mk_signature($z);
			$drug_class_name{$zsig} = $z;
			$$visited{$zsig} = 1;
		}
	}
	
	return sort grep { defined } map { ( $drug_preferred_name{$_} // $drug_class_name{$_} ) } keys %visited_stack;
}


sub get_preferred_drug_name {
# 	print STDERR join(" ", caller() )."\n";
	&ON_DEMAND_INIT;
	if ( exists $drug_preferred_name{ mk_signature($_[0]) } ) {
		return ( $drug_preferred_name{ mk_signature($_[0]) } );
	}
	$drugs_not_found{ $_[0] } ++ if defined $_[0] ;
	return $_[0];
}

my %normalised_treatment_class_name_cache;
sub get_normalised_treatment_class_name { 
	my $x = shift;
	
	return undef if ! defined $x;

	if ( exists $normalised_treatment_class_name_cache{$x} ) {
		if ( ! defined $normalised_treatment_class_name_cache{$x} ) {
			warn "Therapy.pm: normalised_treatment_class_name_cache\{$x\} is undefined." ;
			return undef;
		} 
		return $normalised_treatment_class_name_cache{$x} ;
	}
	
	&ON_DEMAND_INIT;
	
	if ( $x =~ /(?: +\+ +|\|)/ ) {
		my @parts = split /(?: +\+ +|\|)/, $x;
		return ($normalised_treatment_class_name_cache{$x} = join(" + ", (sort ( map { get_normalised_treatment_class_name($_) } @parts ) ) ) );
	}
	
	return ($normalised_treatment_class_name_cache{$x} = ( $drug_class_name{ mk_signature($x) } // $x ) );
}

sub get_drug_synonyms {
	&ON_DEMAND_INIT;
	return sort keys %{ $drug_synonyms{ mk_signature( get_normalised_treatment_name($_[0]) ) } };
}


# print "<[$arg]\t".join("|", (keys %{ $drug_synonyms{ $drug_preferred_name{ mk_signature( get_preferred_drug_name($arg) ) } } } ) )."\n";
# 
our $sep_regex = "[[:space:]\\-]?";
# our $sep_regex_qe = "\\\\Q$sep_regex\\\E";
# 
# sub mk_regex_search_tree(\@;$\%) {
# 	my @x = @{ $_[0] };
# 	my $func = $_[1];
# 	my $previndexref = $_[2];
# 	my $regex; 
# 	my %index;
# 	
# 	if (defined $previndexref) {
# 		for my $k (keys %{ $previndexref }) {
# 			next if ! length $k;
# 			my ($a,$b) = ( $k =~ /^(.)(.*)/ );
# 			print STDERR "$k -- $a-$b ".$$previndexref{$k}."\n";
# 			die if ! defined $a;
# 			$index{$a}{$b} = $$previndexref{$k};
# 		}
# 	} else {
# 		for (@x) {
# 			my ($a,$b) = /^($sep_regex_qe|.)(.*)/;
# 			$index{$a}{$b} = $_;
# 		}
# 	}
# 	
# 	my @elems;
# 	
# 	for my $first (keys %index) {
# 		my @k = sort { length $b <=> length $a } ( keys %{ $index{$first} } ) ;
# 		my $z = '';
# # 		print join("-", @k)."\n";
# 		if ( (scalar(@k) > 2) and (length($k[0]) > 5) ) {
# 			my $zero_length = undef;
# 			if ( grep { ! length } @k ) {
# 				$zero_length = ( $func ? $func->('', $index{$first}{''}) : $index{$first}{''} );
# 				@k = grep { length } @k ;
# 			}
# 			$z = "(?:$first(?:".&mk_regex_search_tree(\@x, undef, \%{ $index{$first} } ).(defined $zero_length ? "|$zero_length" : '')."))";
# 		} else {
# 			$z = "(?:$first(?:". join('|', ( map { $func ? $func->($_, $index{$first}{$_}) : $_ } @k ) ) ."))";
# 		}
# 		push @elems, $z;
# 	}
# 	
# 	return  join("|", @elems);
# }

sub get_drug_regex {
	&ON_DEMAND_INIT;
	my @elems ;
	
	for my $arg (@_) {
		push @elems, map { s/\.(alpha|beta|gamma|delta)\./$1/i; $_ } grep { ! /\d+,\d+-|'|[\[\]]|^\.|, / } keys %{ $drug_synonyms{ mk_signature( get_preferred_drug_name($arg) ) } };
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

# 	my $regex = mk_regex_search_tree(@elems);
	
	my $regex = join("|", sort { length $b <=> length $a } @elems );

	return $regex ;
}



sub get_drug_regex_M {
	&ON_DEMAND_INIT;
	my @elems ;
	for my $arg (@_) {
		push @elems, map { s/\.(alpha|beta|gamma|delta)\./$1/i; $_ } grep { ! /\d+,\d+-|'|[\[\]]|^\.|, / } keys %{ $drug_synonyms{ mk_signature( get_preferred_drug_name($arg) ) } };
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
		$pref_name{$_} = get_preferred_drug_name($_) for @m ;
		if (scalar @a >= 2 ) {
			my ($s) = sort { length $b <=> length $a } grep { /[\- ]/ } @a;
			if ( defined ($s) ) {
				$s =~ s/[ \-]/$sep_regex/g;
				my $d = grep { /^$s$/i } @a ;
				if ( $d == scalar(@a) ) {
					@m = ($s);
					$pref_name{$s} = get_preferred_drug_name($a[0]);
				} 
			}
		} 
		push @elems, @m;
	}

# 	my $regex = mk_regex_search_tree(@elems, sub { my $x = $_[0]; my $y = $_[1]; return "$x(?{push \@M,'".$pref_name{$y}."'})" } );

	my $regex = join("|", map { "$_(?{push \@M,'".$pref_name{$_}."'})" } sort { length $b <=> length $a } @elems);
	
	return $regex ;
}


sub get_drug_classes { 
	&ON_DEMAND_INIT;
	my @a ;
	return @a if ! exists $drug_class{ mk_signature($_[0]) };
	@a = sort map { $drug_class_name{$_} } keys %{ $drug_class{ mk_signature($_[0]) } };
	@a = uniq(@a);
# 	print "\e[1;33m".join("\t", @a)."\e[0m\n";
	return @a ;
}




#################################################################################

sub get_combo_classes { # from drugs
	&ON_DEMAND_INIT;
	my @combo_parts = split /\s*(?:\+|\|)\s*/, $_[0];
	
	my @retval;
	my @drug_classes;
	
	for my $c (@combo_parts) {
		my @a = get_drug_classes($c);
		push @drug_classes, \@a;
	}
	
	my $combo = CombinationCounter->new(\@drug_classes);
	while ( ! $combo->overflow() ) {
		my @combo = $combo->get();
		push @retval, join(" + ", ( sort ( uniq( map { get_normalised_treatment_class_name($_) } @combo) ) ) );
		$combo->incr();
	}
	return @retval;
}


sub is_a_treatment_class {
	&ON_DEMAND_INIT;
	my @combo_parts = split /\s*(?:\+|\|)\s*/, $_[0];
	return 0 if ! scalar @combo_parts ;
	for my $c (@combo_parts) {
		return 0 if ! is_a_drug_class($c);
	}
	return 1;
}


sub get_all_drugs { 
	&ON_DEMAND_INIT;
	return sort map { $drug_preferred_name{$_} } keys %drug_class;
}

sub get_all_drug_classes { 
	&ON_DEMAND_INIT;
	return sort map { $drug_class_name{$_} } keys %drug_class_name;
}

sub get_all_known_drugs_for_class { # drug class 
	&ON_DEMAND_INIT;
	my $dc_sig = mk_signature($_[0]) ;
	
	my %retval ; 
	for my $d_sig ( keys %drug_class ) {
		next if ! exists $drug_class{$d_sig}{ $dc_sig };
		$retval{ $drug_preferred_name{$d_sig} } = 1;
# 		print "$dc_sig\n";
	}
	if ( exists $drug_class_is_a{$dc_sig} ) {
		for my $dc ( keys %{ $drug_class_is_a{$dc_sig} } ) {
# 			print STDERR "$dc_sig\t$dc\n";
			my @retval_subclass = get_all_known_drugs_for_class($dc);
			$retval{$_} = 1 for @retval_subclass ;
		}
	}
	return sort keys %retval ; 
}

sub get_all_drug_class_pairs {
	&ON_DEMAND_INIT;
	my @pairs;
	for my $d_sig ( keys %drug_class ) {
		for my $dc_sig ( keys %{ $drug_class{$d_sig} } ) {
			push @pairs, join("\t", get_preferred_drug_name($d_sig), get_normalised_treatment_class_name($dc_sig) );
		}
	}
	return @pairs;
}

sub get_all_drug_class_hierarchy {
	&ON_DEMAND_INIT;
	return @hierarchy_pairs;
}

sub load_drug_database {
	my $database_file = shift;
	my @lines = file($database_file);
	
	for my $line ( @lines ) {
		chomp $line;
		$line =~ s/#.*// ; # remove comments;
		next if $line =~ /(?:^|\b)drug_class(?:\b|$)|^\s*$/ ;  # get rid of header
		
		my ($drug_str, $drug_classes) = split /\t/, $line;
		my @drug_parts = map { trim($_) } split /\|/, $drug_str;
		my $primary_drug_name = $drug_parts[0];

		if ( $drug_classes =~ / \+ / ) { # this is a drug combination regimen
			my $drug_combinations = $drug_classes;
			$drug_combo_regimens{ mk_signature($_) } = trim($drug_combinations) for @drug_parts;
			next;
		}
		
		for my $drug (@drug_parts) {
			if ( ! exists $drug_preferred_name{ mk_signature($drug) } ) {
				$drug_preferred_name{ mk_signature($drug) } = $primary_drug_name ;
				$drug_name{ mk_signature($drug) } = $primary_drug_name ;
			}
			
			$drug_synonyms{ mk_signature($primary_drug_name) }{ $drug } ++;
			
			warn "\e[1;31mDATABASE WARNING: $database_file: Drug ``$drug_str'' has no drug classes.\e[0m\n" if ! defined $drug_classes ;
			my @dc_parts = map { trim($_) } split/\s*;\s*/, $drug_classes;
			for my $dcp ( @dc_parts ) {
				$drug_class{ mk_signature($drug) }{ mk_signature($dcp) } = 1;
				$drug_class_name{ mk_signature($dcp) } = $dcp ;
			}
		}
	}
}

sub load_drug_databases {
	return %drug_class if scalar %drug_class;
	
	my @drug_database_files = @_;
	
	for my $db (@drug_database_files) {
		if ( ! -f $db ) {
			warn "\e[1;31mDATABASE WARNING: $db not found.\e[0m\n" ;
			next;
		}
		load_drug_database($db);
	}
}

sub load_drug_class_hirerchy {
	my $drug_class_hirarechy_file = shift;
	for ( file( $drug_class_hirarechy_file ) ) {
		chomp;
		my @parts = map { trim($_) } split /\t/, $_;
		for my $p (@parts) {
			next if exists $drug_class{ mk_signature($p) }; # this is a drug, not a drug class
			$drug_class_name{ mk_signature($p) } = $p;
		}
		for my $i (0 .. $#parts) {
			last if $i+1 > $#parts;
# 			print STDERR " mk_signature($parts[$i]) =>  mk_signature($parts[$i+1]) \n";
# 			print STDERR "\e[1;31mTherapy.pm: Error in drug class hierarchy: $parts[$i] == $parts[$i+1]\e[0m" if $parts[$i] eq $parts[$i+1];
			$drug_class_is_a{ mk_signature($parts[$i]) }{ mk_signature($parts[$i+1]) } = 1;
			$drug_class_has_parent{ mk_signature($parts[$i+1]) }{ mk_signature($parts[$i]) } = 1;
			push @hierarchy_pairs, join("\t", $parts[$i], $parts[$i+1]);
		}
	}
	
	$drug_regex = get_drug_regex( get_all_drugs() );
	
	return %drug_class ;
}


sub dc_lev { # normalised levenshtein
	my ($x, $focus) = ( $_[0], $_[1] );
# 	print join(" ", ($x, $focus) )."\n";
	return levenshtein($focus, $x) / ( length($x) + length($focus) + 1);
}

sub lcss { # longest common substring 
	my @strings = sort { length $a <=> length $b } @_;
	my ($x, $y) = @strings ;
	my $xl = length $x;
	my $lcss = '';
	for ( my $i=0; $i<$xl; ++$i ) {
		for ( my $j=length($lcss)+1; ($i+$j)<=$xl; ++$j ) {
			my $xss = substr $x, $i, $j;
			last if ( index($y, $xss) == -1 );
			$lcss = $xss ;
		}
	}
	return $lcss;
}

our %normalised_treatment_name_cache = ();

sub get_normalised_treatment_name {
	my $treatment = shift;
	return $normalised_treatment_name_cache{$treatment}  if exists $normalised_treatment_name_cache{$treatment} ;
	&ON_DEMAND_INIT;
	my $sig_treatment = mk_signature($treatment) ;
	$treatment = $drug_combo_regimens{ $sig_treatment } if exists $drug_combo_regimens{ $sig_treatment };
	my @combo_parts = map { trim($_) } split /\s*(?:\||\+)\s*/, $treatment;
	@combo_parts = sort map { get_preferred_drug_name($_) } @combo_parts ;
	return ( $normalised_treatment_name_cache{$treatment} = join(" + ", @combo_parts ) );
}

sub get_treatment_class {
	&ON_DEMAND_INIT;
	
	my $treatment = shift;
	
	if ( is_a_treatment_class($treatment) ) {
# 		print STDERR "\e[0;31mTherapy.pm: $treatment is already a treatment class.\e[0m\n";
		return get_normalised_treatment_class_name($treatment);
	}
	
	my @focus = grep { defined } @_; # a list of strings of biomarker names on a particular class a drug has mutliple targets
	warn join("\t", caller)."\n" if ! defined $treatment;
	
	my @combo_parts = map { trim($_) } split /\s*\+\s*/, $treatment;
	
	@combo_parts = map {
		my $dc = $_;
		my $dc_sig = mk_signature( $dc );
		my @retval = ($dc);
		if ( exists $drug_combo_regimens{$dc_sig} ) {
			my @drugs = split /\s+\+\s*/, $drug_combo_regimens{ $dc_sig };
			@retval = @drugs;
			}
		@retval
		} @combo_parts ;

	my @combo_parts_class ;
	for my $part (@combo_parts) {
		my @dcs = get_drug_classes($part); 

		my %dcs_min;
		for my $dc (@dcs) {
			for my $focus (@focus) {
				my $lcssl = length( lcss( $dc, $focus ) );
# 				printf "$dc\t$focus\t$lcssl\n";
				if ( ! exists $dcs_min{$dc} ) {
					$dcs_min{$dc} = $lcssl ;
				} else {
					$dcs_min{$dc} = max( $dcs_min{$dc}, $lcssl ) ;
				}
			}
		}
# 				@dcs = sort { 
# 				length( lcss( $b, $focus ) ) <=> length ( lcss( $a, $focus ) ) ||
# 				length $b <=> length $a
# 				} @dcs ;
				
# 			if ( grep { /$focus/ } @dcs ) {
# # 				@dcs = sort { dc_lev( $a, $focus ) <=> dc_lev( $b, $focus ) } grep { /$focus/ } @dcs ;
# 				@dcs = sort { length( lcss( $b, $focus ) ) <=> length ( lcss( $a, $focus ) ) } grep { /$focus/ } @dcs ;
# 			} else {
# 				@dcs = sort { length( lcss( $b, $focus ) ) <=> length ( lcss( $a, $focus ) ) } @dcs ;
# # 				@dcs = sort { length $b <=> length $a } @dcs ;
# 			}
#			@dcs = map { dc_lev($_, $focus)." $_" } sort { dc_lev( $a, $focus ) <=> dc_lev( $b, $focus ) } @dcs ;
		
		if ( scalar %dcs_min ) {
			@dcs = sort { ($dcs_min{$b} <=> $dcs_min{$a}) || (length($a) - length($b) )} keys %dcs_min ;
		}
		
		push @combo_parts_class, ( scalar(@dcs) ? $dcs[0] : 'UNDEF' ); # $dcs[0]
	}
	return join(" + ", ( sort map { get_normalised_treatment_class_name($_) } @combo_parts_class) );
}

sub is_combination {
	return ( $_[0] =~ /(?:\s+\+\s+|\|)/ );
}

sub register_known_combination {
	if ( is_combination($_[0]) ) {
		$known_drug_combinations{ get_normalised_treatment_name($_[0]) } = 1;
	}
}

sub get_all_known_combinations {
	return sort keys %known_drug_combinations;
}


my %treatment_contains_drug_cache;

sub treatment_contains_drug {
	my $treatment = shift;
	my $drug = shift;
	return $treatment_contains_drug_cache{$treatment}{$drug} if exists $treatment_contains_drug_cache{$treatment} and exists $treatment_contains_drug_cache{$treatment}{$drug} ;
	
	&ON_DEMAND_INIT;
	
	my @combo_parts = map { trim($_) } split /\s*(?:\+|\|)\s*/, $treatment;
# 	print join("\t", 'treatment_contains_drug', $treatment, $drug)."\n";
	for my $p (@combo_parts) {
		if ( get_preferred_drug_name($p) eq get_preferred_drug_name($drug) ) {
			$treatment_contains_drug_cache{$treatment}{$drug} = 1;
			return 1;
		}
	}
	$treatment_contains_drug_cache{$treatment}{$drug} = 0;
	return 0;
}

sub get_all_drug_class_offspring($;$) {
	my $dc = shift;
	my $visited = shift;
	
	my %visited;
	$visited = \%visited if ! defined $visited;
	
	my $dcsig = mk_signature($dc) ;
	
	return if exists $visited{$dcsig};
	$visited{$dcsig} = 1;
	
	my @retval;
	push @retval, $drug_class_name{$dcsig};
	
# 	print $dcsig."\t".$drug_class_name{$dcsig}."\n";
	if ( exists $drug_class_is_a{$dcsig} ) {
		my @children = sort keys %{ $drug_class_is_a{ $dcsig } };
		for my $c (@children) {
			push @retval, &get_all_drug_class_offspring($c, $visited);
		}
	}
	
	return @retval;
}


sub get_all_matched_offspring_treatment_classes($) {
	&ON_DEMAND_INIT;
	my $treatment_class = $_[0];
# 	my @dc_combo_parts = map { exists $drug_combo_regimens{ mk_signature($_)} ? ( split /\s*(?:\+|\|)\s*/, $drug_combo_regimens{ mk_signature($_) } ) : trim($_) } split /\s*(?:\+|\|)\s*/, $treatment_class;
	my @dc_combo_parts = map { trim($_) } split /\s*(?:\+|\|)\s*/, $treatment_class;
	my @dc_combo_parts_expanded;
# 	print STDERR ">> ".(join(" + ", @dc_combo_parts ))."\n";
	for my $dc ( @dc_combo_parts ) {
		my @offsprings = grep { defined } get_all_drug_class_offspring($dc);
		push @dc_combo_parts_expanded, \@offsprings;
# 		print STDERR "$dc: ".(join("; ", @offsprings))."\n";
	}
	
	my $combo = CombinationCounter->new(\@dc_combo_parts_expanded);
	my @retval;
	while ( ! $combo->overflow() ) {
		my @combo = $combo->get();
# 		print STDERR join(" ++ ", map { defined ? $_ : "\e[1;36mUNDEFINED\e[0m" } @combo)."\n";
		push @retval, join(" + ", ( sort map { get_normalised_treatment_class_name($_) // 'UNDEF' } @combo) );
		$combo->incr();
	}
	return sort @retval;
	
}


sub ON_DEMAND_INIT {
# 	print "!";
	return if $f_initiailised ;
	$f_initiailised = 1;
	
	for my $srcfile ( POTTRConfig::get_paths('data', 'therapy-database-file') ) {
		load_drug_databases $srcfile ;
	}

	for my $srcfile ( POTTRConfig::get_paths('data', 'drug-class-hierarchy-file') ) {
		load_drug_class_hirerchy $srcfile ;
	}
}


##########################################################################################################################
#
# Examples - Loading drug database
#
# print get_treatment_class('Vemurafenib')."\n";
# print get_treatment_class('Vemurafenib+Cobimetinib')."\n";
# print get_treatment_class('dd+Cobimetinib')."\n";
# 
# print get_preferred_drug_name('PLX4032')."\n";
# print get_preferred_drug_name('veMURAFenib')."\n";
# 
# print map { "Synonyms: $_\n" } get_drug_synonyms('veMURAFenib');
# print map { "Class: $_\n" } get_drug_classes("Regorafenib");
# print map { "All drug in class: $_\n" } get_all_known_drugs_for_class('ALK_inhibitor');
# print map { "All drug in class: $_\n" } get_all_known_drugs_for_class('MET_inhibitor');
# 
# print is_a('ALK_inhibitor', 'ALK_inhibitor,first_generation' )."\n";
# print is_a('Selumetinib', 'MEK_inhibitor')."\n";
# print is_a('Placebo', 'Placebo' )."\n";
#
#
# # for my $drug ( get_all_drugs() ) {
# # 	print "$drug:\n";
# # 	print "\t$_\n" for get_drug_classes($drug);
# # }
# 
# for my $drug ( get_all_known_drugs_for_class('immune_checkpoint_blockade') ) { # immune checkpoint blockade,PD-1 targeting
# 	print "$drug:\n";
# 	print "\t$_\n" for get_drug_classes($drug);
# }
# 
# print is_a_drug_class("MEK_inhibitor")."\n";
# print is_a_drug_class("Vemurafenib")."\n";
# 
# print sort map { "$_\n" } get_all_parents("Pembrolizumab");
#
# print sort map { "$_\n" } get_all_drug_class_hierarchy();
# 
# print sort map { "$_\n" } get_all_drug_class_pairs();
#
# print get_drug_regex("Vemurafenib", "Cobimetinib");

# print get_drug_regex("Vemurafenib");
# print get_drug_regex( get_all_drugs() );
1;
