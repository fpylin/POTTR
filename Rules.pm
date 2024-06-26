#!/usr/bin/perl
##############################################################################
#
# Rules.pm - a rule-based, forward-chaining engine with inference tracking
# 
# This file is part of Oncology Treatment and Trials Recommender (OTTR) and 
# Precision OTTR (POTTR).  The history of facts derived  from inference are 
# automatically tracked to permit retrospective auditing and analysis.
# 
# SYNOPSIS:
# 
# Creating a new ruleset:
# 
# my $ruleset = Rules->new( 
#  { 'rules' => \@rules, 
#    'debug_level' => 0, 
#    'debug_msg' => 'ansi' 
#  } 
# ) ;
# 
# my $new_rule = Rules::mkrule( ['A', 'B'], ['C'] );  # A; B => C
# 
# my @rules = ($new_rule, ...);
# 
# $ruleset ->load(@rules);
#
# my $facts = Facts -> new(@facts);
# my %facts_with_tags = $ruleset ->run( $facts ) -> get_fact_tag_hash();
# 
# $ruleset ->debug_print(); # prints debug messages to STDIN;
# $ruleset ->to_string();
# $ruleset ->print();
# 
# Two types of rules:
# 
# 1. STATIC RULES are described using the following format:
#   A => B              # if A then B
#   A; B => C           # if A AND B then C
#   A; NOT B => C       # if A AND NOT B then C
#   A => B [TAG1]       # if A then B, with B now tagged with TAG1
#   A; TRACK TAG1 => B  # same as above, allows inference tracking
# 
# 2. DYNAMIC RULES are evaluated by calling a perl subroutine, such that:
#
# $rules->define_dyn_rule('rule-name', \&rule_subroutine)
#
# sub rule_subroutine {
# 	my $triggering_fact = $_[0];
# 	my $facts_hash_ref  = $_[1];
#   # ... 
# 	return ('New Fact 1');
# }
# 
# 3. The Facts class may contain variables and function that are not part 
#    of the production rules
##############################################################################
#
# VERSION HISTORY 
# 
# The inference engine is derived from ESCORT.pl - an expert system prototype 
# (Expert System for evaluating Clinical Outcomes and Response to Treatment) 
# that  uses  the  Cancer Care: Treatment Outcome (CCTO) ontology for inferring 
# high-level concepts about treatment outcomes of cancer cases from a list  
# of cancer-related clinical events. 10 December 2017
# 
# Revision 0.1   - 22 May 2019
# Revision 0.8   - 27 Nov 2019
# Revision 0.9   - 18 Dec 2019
# Revision 0.91  - 04 May 2019
# Revision 0.92  - 13 Jun 2019
# Revision 0.95  -  2 Jan 2022
# 
# Copyright 2019-2021, Frank Lin & Kinghorn Centre for Clinical Genomics, 
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

{ # class / package for managing facts
package Facts;

use Storable;

sub uniq { my %a; $a{$_} = 1 for(@_) ; return sort keys %a; }

sub new {
    my $class = shift;
    my @facts = @_;
    
    my %a;
    my $self = { 'facts' => \%a };

    bless($self, $class);
    
	for my $s (@facts) {
		$self->assert_string($s);
	}
	
	my %hash_data;
	my %hash_func;
	
	$self->{'data'} = \%hash_data;
	$self->{'func'} = \%hash_func ;
	
    return $self;
}

my %escchar_map = (
# 	"\n" => 'LF',
	';' => 'SEMICOLON',
	'[' => 'LSQ',
	'#' => 'HASH',
	']' => 'RSQ'
	);
	
my %escchar_map_r = map { $escchar_map{$_} => $_ } keys %escchar_map ;
my $escchar_map_r_regex = join('|', sort { length $b <=> length $a } keys %escchar_map_r );

my %__escape_cache;

sub escape {
	my $x0 = my $x = shift;
	return $__escape_cache{$x} if exists $__escape_cache{$x} ;
	$x =~ s/([;#\[\]])/"--".$escchar_map{$1}."--"/sge;
	$__escape_cache{$x0} = $x ;
	return $x;
}

sub unescape {
	my $x = shift;
	$x =~ s/--($escchar_map_r_regex)--/$escchar_map_r{$1}/sge;
	return $x;
}

our %mk_fact_str_cache ;

sub mk_fact_str($@) {
	my $fact = shift;
	my @tags = @_;
	@tags = grep { defined } @tags ;
	return $fact if ! scalar @tags;
	my $z = join("\0",$fact, @tags) ;
	return $mk_fact_str_cache{ $z } if ( exists $mk_fact_str_cache{ $z } );
	@tags = uniq( @tags );
	return (  $mk_fact_str_cache{ $z } = $fact." [".join('; ', ( map { escape($_) } @tags))."]"  );
}

our %string_to_fact_and_tags_cache;

sub string_to_fact_and_tags($) {
	my $input = shift; # string
	return @{ Storable::thaw( $string_to_fact_and_tags_cache{$input} ) } if exists $string_to_fact_and_tags_cache{$input};
	my ($fact, $tags) = ( $input =~ /^(.+)\s*\[\s*(.+)\s*\]\s*$/ );
	$fact = $input if ! defined $fact;
	$fact =~ s/^\s*//; $fact =~ s/\s*$//;
	my @tags ;
	@tags = split /\s*;\s*/, $tags if defined $tags ; # map { unescape($_) } 
	my @a = ($fact, @tags);
	$string_to_fact_and_tags_cache{$input} = Storable::freeze(\@a);
	return @a;
}

sub string_to_tags_hashed($) {
	my $input = shift; # string
	my ($tags) = ( $input =~ /^.+\s*\[\s*(.+)\s*\]\s*$/ );
	my @tags ;
	@tags = split /\s*;\s*/, $tags if defined $tags ; # map { unescape($_) }
	my %tags ;
	for my $t (@tags) { 
		my ($h, $v) = split /:/, $t, 2; 
		my $vv = unescape($v // ''); 
		if ( ! exists $tags{$h} ) {
			$tags{$h} = $vv ;
		} else {
			$tags{$h} .= "; $vv" ;
		}
	} 
	return %tags;
}

sub assert_tags { 
	my $self = shift;
	my $fact = shift;
	my @tags = @_;
# 	print ">>   $fact\t".join("; ", keys %{ $self->{'facts'}{$fact} } )."\n";
	$self->{'facts'}{$fact}{$_} = 1 for @tags ;
# 	print ">>>> $fact\t".join("; ", keys %{ $self->{'facts'}{$fact} } )."\n";
}

sub assert { # returns if assertion is a new fact
	my $self = shift;
	my $fact = shift;
	my @tags = @_;
	my $retval = 1;

	if ( exists $self->{'facts'}{$fact} ) {
		$retval = 0 ;
	} else {
		my %a ;
		$self->{'facts'}{$fact} = \%a;
	}
	$self->assert_tags($fact, @tags);
	
	return $retval;
}

sub retract { # retract a fact
	my $self = shift;
	my $fact = shift;
	my $retval = 0;

	if ( exists $self->{'facts'}{$fact} ) {
		delete $self->{'facts'}{$fact};
		$retval = 1 ;
	} 
	
	return $retval;
}

sub assert_string { # returns if assertion is a new fact
	my $self = shift;
	my $input = shift;
	my ($f, @tags) = string_to_fact_and_tags($input);
	return $self->assert($f, @tags);
}

sub has_fact { 
	my $self = shift;
	my $f = shift;
	return exists( $self->{'facts'}{$f} );
}

sub get { # return the hash of facts/tags
	my $self = shift;
	return $self->{'facts'} ;
}

sub get_facts_list { # return a list of facts
	my $self = shift;
	return sort keys %{ $self->{'facts'} } ;
}

sub get_tags {
	my $self = shift;
	my $f = shift;
	if ( exists $self->{'facts'}{$f} ) {
		return keys %{ $self->{'facts'}{$f} } ;
	} else {
		my @null;
		return @null ;
	}
}

sub untag {
	my $self = shift;
	my $f = shift;
	my @tags_to_remove = @_;
	if ( exists $self->{'facts'}{$f} ) {
		for my $t (@tags_to_remove) {
			delete $self->{'facts'}{$f}{$t} if exists $self->{'facts'}{$f}{$t};
		}
	}
}

sub get_fact_tag_hash { # return a get_fact_tag_hash 
	my $self = shift;
	my @facts = $self->get_facts_list();
	return map { $_ => join("; ", sort( $self->get_tags($_) ) ) } @facts  ;
}

sub get_fact_tag_strings { # return a list of facts
	my $self = shift;
	my @facts = $self->get_facts_list();
	return map { 
		my $f = $_;
		my @tags = $self->get_tags($f);
		mk_fact_str($f, @tags)
		} @facts ;
	return sort keys %{ $self->{'facts'} } ;
}

sub to_strings {
	my $self = shift;
	return sort keys %{ $self->{'facts'} } ;
}

sub clean_tags { # SUBROUTINE TO DELETE
	my $self = shift;
	my @facts = $self->get_facts_list();
	
	for my $f (@facts) {
		my %to_delete;
		my @trk = sort { length($a) <=> length($b) } ( $self->get_tags($f) );
		TRKI: for my $i (0..$#trk) {
			for my $j ( ($i+1) ..$#trk) {
				if ( index($trk[$j], $trk[$i]) != -1 ) {
# 					print "!A $trk[$i]!\n";
					$to_delete{ $trk[$i] } = 1;
					next TRKI;
				}
			}
		}
		@trk = grep { ! exists $to_delete{$_} } @trk;
		my %trk = map { $_ => 1 } @trk;
		$self->{'facts'}{$f} = \%trk;
	}
}

# the Facts class also include the capability some global variables

sub setvar($$) { # set "global" variables
	my $self = shift;
	my $f = shift;
	my $v = shift;
# 	print STDERR "Rules.pm: Setting $f = $v\n";
	$self->{'_data_'}{$f} = $v;
}

sub deffunc($\&) { # set "global" functions
	my $self = shift;
	my $f = shift;
	my $func_ref = shift;
# 	print STDERR "$f\t$func_ref\n";
	$self->{'_func_'}{$f} = $func_ref;
}

sub getvar($) {
	my $self = shift;
	my $f = shift;
# 	print STDERR "$f\t$self->{'data'}{$f}\n";
	return undef if ! exists $self->{'_data_'}{$f};
	return $self->{'_data_'}{$f};
}

sub func($@) {
	my $self = shift;
	my $f = shift;
	my @args = @_;
	if ( ! exists $self->{'_func_'}{$f} ) {
		warn "Rules.pm: function $f is not defined";
		return undef ;
	}
	return $self->{'_func_'}{$f}->($self, @args);
}

1;
};

######################################################################################
{
package Rules;

use Storable;

our @EXPORT;
our @EXPORT_OK;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 0.8 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw( &mkrule &facts_to_graphviz );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}

sub mkrule { 
	my $lhs = shift;
	my $rhs = shift;
	my $attrib = shift;
	my $rule = ( defined $attrib ?  "[".join('; ', ( grep { defined and length } @{ $attrib } ) )."] " : '' ).
		join(' => ', 
			join('; ', ( grep { defined and length } @{ $lhs } ) ), 
			join('; ', ( grep { defined and length } @{ $rhs } ) ) 
		); 
	return $rule;
}

sub new {
    my $class = shift;
    my %params ;
    %params = %{ $_[0] } if defined $_[0] ;
    
    my @a;
    my $self = { 
		'rules' => \@a, 
		'debug_level' => 0, 
		'debug_msg' => 'ansi',
		'run_rules_only_once' => 1, # if specified, run rules more than once.
	};
    
    bless($self, $class);
    
    ${$self}{$_} = $params{$_} for (keys %params) ;
    
    $self->{'f_needs_reindex'} = 0;
    
    return $self;
}


sub debug_print {
	my ($self, $color, @msgs) = @_;
	my $debug_output = '';
	for my $msg (@msgs) {
		for ( $self->{'debug_msg'} ) {
			/ansi/ and do { $debug_output .= "\e[1;${color}m".$msg."\e[0m"."\n"; }; 
			/html/ and do { $debug_output .= "<span class=c${color}m>".$msg."</span><br/>\n"; };
		}
	}
	print $debug_output if $self->{'debug_level'} > 0;
	$self->{'debug_output'} .= $debug_output ;
}

sub index_rules { # Rules::index_rules - index the rules by lhs patterns
	my $self = shift;
	for my $i ( 0 .. $#{ $self->{'rules'} } ) {
		for my $j ( 0 .. $#{ $self->{'rules'}[$i]{'lhs'} } ) {
			my $lhse_neg = $self->{'rules'}[$i]{'lhs_neg'}[$j] ;
			my $lhse = $self->{'rules'}[$i]{'lhs'}[$j];
			next if $lhse_neg != 0 ;
			$self->{'rule_index'}{$lhse}{$i} = 1;
# 			print "\@ $lhse <= [$i] ".&rule_to_string( $self->{'rules'}[$i] )." \n";
		}
	}
    $self->{'f_needs_reindex'} = 0;
}


sub retrieve_rules { # Rules::retrieve_rules - index the rules by lhs patterns
	my ($self, @facts) = @_;
	my %ruleset ;
	for my $fact (@facts) {
		next if ! exists $self->{'rule_index'}{$fact};
		for my $i ( keys %{ $self->{'rule_index'}{$fact} } ) {
			$ruleset{$i} ++; 
		}
	}
	return map { $self->{'rules'}[$_] } keys %ruleset ;
}


sub retrieve_rules_by_prio { # Rules::retrieve_rules - index the rules by lhs patterns, returns $hash{prio_score} = \@rules
	my ($self, @facts) = @_;
	my %ruleset ;
	for my $fact (@facts) {
# 		print "\t retrieve_rules: $fact\n";
		next if ! exists $self->{'rule_index'}{$fact};
		for my $i ( keys %{ $self->{'rule_index'}{$fact} } ) {
# 			print "| $fact => [$i] ".&rule_to_string( $self->{'rules'}[$i] )." \n";
			$ruleset{$i} ++; 
		}
	}
	my %retval;
	
	for my $rid ( keys %ruleset ) {
		my $prio = $self->{'rules'}[$rid]{'prio'} ;
		push @{ $retval{$prio} }, $self->{'rules'}[$rid];
	}
	
	for my $prio ( keys %retval ) {
		my @ary = @{ $retval{$prio} };
		@ary = sort { $$a{'order'} <=> $$b{'order'} } @ary;
		$retval{$prio} = \@ary;
	}
	
	return %retval;
}


sub rule_to_string {
	my $rule = $_[0];
	my $output = 
		($rule->{'prio'} ? "[prio:$rule->{'prio'}] " : '').
		join("; ", map { ($rule->{'lhs_assert'}[$_]?'*':'').($rule->{'lhs_neg'}[$_]?'NOT ':'').$rule->{'lhs'}[$_] } (0 .. $#{$rule->{'lhs'}} ) ).
		" => ".
		join("; ", map { ($rule->{'rhs_neg'}[$_]?'NOT ':'').$rule->{'rhs'}[$_] } (0 .. $#{$rule->{'rhs'}} ) )
		;
	return $output;
}

sub to_string { 
	my ($self, @ruleset) = @_ ;
	@ruleset = @{ $self->{'rules'} } if ! scalar (@ruleset) ;
	my $output;
	
	for my $i ( 0..$#ruleset ) {
		$output .= join("\t", "r-$i", rule_to_string($ruleset[$i]))."\n";
	}
	return $output // '' ;
}
	

sub print {
	my $self = shift;
	print $self->to_string(@_);
}

sub clean_fact {
    my $fact = shift;
    my $a = $fact;
    $a =~ s/\s*\[.+\]//;
    return $a;
}


sub fire_rule {
	my $self = shift; 
	my $facts = shift;
	my $rule = shift; 
	my $constraints = shift; # optional reference to additional tracking tags ;
	my $n_new = 0;
	my @rhs = @{ $$rule{'rhs'} };
	my @rhs_neg = @{ $$rule{'rhs_neg'} };
	return 0 if ! scalar (@rhs);

	for my $i (0 .. $#rhs) {
		my $e = $rhs[$i];
		my @tags;
		($e, @tags) = Facts::string_to_fact_and_tags( $e );
		my $n_tags = scalar @tags;
		@tags = grep { $_ ne '__untrack__' } @tags; # if defined __untrack__ then no history should be tracked.
		my $f_untrack = ( $n_tags != scalar(@tags) ) ;
		
		if ( $rhs_neg[$i] ) { # NOT in RHS indicates that a fact should be retracted
			if ( $facts->has_fact($e) ) { 
				$n_new += $facts->retract($e); 
			} 
			next;
		} else { # assign the fact specified in the RHS
			if ( ! $facts->has_fact($e) ) { # if fact already exists but no new tag defined, then simply do nothing.
				$n_new += $facts->assert($e, @tags); # if a fact does not exist, then add a new fact!
			}
		}
			
		my @lhs_tags;
		
		if (! $f_untrack) {
			for my $i ( 0 .. $#{ $$rule{'lhs'} } ) {
				my $lhse = $$rule{'lhs'}[$i];
# 				next if ! $facts->has_fact($lhse); # skip if the fact does not exist CHECK; (or no tags are mentioned in LHS element)
				
				my $lhse_cleaned = clean_fact(($$rule{'lhs_neg'}[$i] ? "NOT " : "").$lhse);
				my @t = $facts->get_tags($lhse);
				push @t, $lhse_cleaned if $lhse_cleaned ne $e;
				push @lhs_tags, @t;
			}
		}
		
		if ( defined $constraints and defined $constraints->{'tags'} )  { # and scalar( @{ $constraints->{'tags'} } )
			my @new_tags = @{ $constraints->{'tags'} } ;
# 			print "! ".scalar( @{ $constraints->{'tags'} } )."\n";
# 			die if scalar( @{ $constraints->{'tags'} } );
# 			print "<>$_<>\n" for @new_tags ;
			push @tags, @new_tags ;
		}
		
		$facts->assert_tags($e, @lhs_tags, @tags);
	}
	return $n_new;
}

sub match_rule_lhs {
	my $self = shift;
	my $facts = shift;
	my $rule = shift;
	my @lhs = @{ $$rule{'lhs'} };
	my @lhs_neg = @{ $$rule{'lhs_neg'} };
	my @lhs_assert = @{ $$rule{'lhs_assert'} };
	my @matched_cons; # matched constraints
	my @unmatched_cons; # unmatched constraints
	my @tracking_tags;
	my @lhs_matched ;
	for my $i (0..$#lhs){ 
		my $lhs_cons = $lhs[$i]; # LHS constraint
		my $f_has_fact ;
		if ( $lhs_cons =~ /^TRACK +(.+)/ ) {
			$f_has_fact = 1;
			push @tracking_tags, $1;
		} else {
			$f_has_fact = $facts->has_fact($lhs_cons) ;
		}
		my $tf = ( ($lhs_neg[$i] != 0) ? # Handling negation: treating it as an XOR problem.
			( $f_has_fact ? 0 : 1 ) :
			( $f_has_fact ? 1 : 0 ) );
		$lhs_matched[$i] = $tf;
		if ($tf) { # constraint i satisfied
			push @matched_cons, $lhs_cons ;
		} else { # constraint i not satisfied
			if ( $lhs_assert[$i] ) { # contingency reasoning by falsely asserting a fact to satisfy a constraint. This fact is tracked to allow retrospective inference
				push @matched_cons, $lhs_cons ;
				my $ttag = "*".($lhs_neg[$i] ? "NOT ": "").$lhs_cons ;
				push @tracking_tags, $ttag ;
				$lhs_matched[$i] = 1;
			} else { # normal behaviour
				push @unmatched_cons, $lhs_cons ;
			}
		}
	}

	my %constraints = (
		'matched' => \@matched_cons,
		'unmatched' => \@unmatched_cons,
		'tags' => \@tracking_tags
	);
	
	my $n_lhs_matched = () = grep { $_ } @lhs_matched ;
	my $n_lhs = scalar(@lhs);
	my $f_matched = ( $n_lhs == $n_lhs_matched );
	return $f_matched, \%constraints ; # ($n_lhs_matched, $n_lhs);
}


sub run { # Rules::run - the main inference engine. 
	my $self = shift;
	my $facts_ref = shift;  # Input: a Facts object
	$self->{'debug_output'} = '';
	
	$$facts_ref->assert('(initial-fact)'); #  push @facts, '(initial-fact)';
	
	$self->index_rules() if $self->{'f_needs_reindex'} ;
	
	my $n_new=0;
	my $iter=0;
	
	my %dyn_rule_facts_matched;
	my @dyn_ruleset = keys %{ $self->{'dyn_rules'} }; # dynamic rule sets;
		
	my %visited_rules;
		
	do	{
		++$iter;
		$n_new=0;

		my @fact_list = $$facts_ref->get_facts_list();
		
		$self->debug_print(32, "Iter $iter: ".scalar(@fact_list). " facts"); # .join(" | ", @fact_list) );
		
		### First run the dynamic rule sets
		my @dyn_rule_unmatched_facts = grep { ! exists $dyn_rule_facts_matched{$_} } @fact_list;
		
		for my $fact (@dyn_rule_unmatched_facts) {
			$dyn_rule_facts_matched{$fact}++;
			for my $dyn_rule (@dyn_ruleset) {
				my @dyn_rule_rhs_to_fire = ( $self->{'dyn_rules'}{$dyn_rule}{'code'}->( $fact, $$facts_ref, \%{ $self->{'cache'} } ) );
				my @dyn_rule_lhs = ($fact);
				my @fake_lhs = ($fact);
				my @fake_rhs_neg ;
				my @fake_lhs_assert = (0);
				my @fake_lhs_neg = (0);
				my %a = ('lhs' => \@fake_lhs, 'lhs_neg' => \@fake_lhs_neg, 'lhs_assert' => \@fake_lhs_assert, 'rhs' => \@dyn_rule_rhs_to_fire, 'rhs_neg' => \@fake_rhs_neg );
				if (scalar @dyn_rule_rhs_to_fire) {
					my $rule = \%a;
				    $self->debug_print(31, "   ".join("  ", $self->{'dyn_rules'}{$dyn_rule}{'name'}, $fact, '=>', join('; ', @dyn_rule_rhs_to_fire)))  ;
					$n_new += $self->fire_rule($$facts_ref, $rule);
# 					print STDERR "\e[1;36m|- ".join("\n", $fact, @{ $a{lhs} }, @{ $a{rhs} })." \e[0m\n" ;
				}
			}
		}
		goto TO_NEXT_CYCLE if $n_new > 0;

		### Now run the static rules clustered by the rule priorities
		my %static_ruleset_by_prio = $self->retrieve_rules_by_prio( $$facts_ref->get_facts_list() );
		
		for my $prio ( sort {$a <=> $b} keys %static_ruleset_by_prio ) {
			for my $rule ( @{ $static_ruleset_by_prio{$prio} } ) { 
				my ( $f_match, $constraints ) = $self->match_rule_lhs($$facts_ref, $rule);
# 				print "<$f_match> $_ <>\n" for @{ $constraints->{tags} } ;
				
				if ( $f_match ) {
					my $match_id = join( "|", $rule->{'order'}, @{ $$constraints{'matched'} } );
# 					$self->debug_print(35, "   visited = ".($visited_rules{$match_id}//0)." mid = ".$match_id )  ;

					if ( ! $self->{'run_rules_only_once'} or ! exists $visited_rules{$match_id} ) {
						my $n_lhs = scalar @{ $$rule{'lhs'} };
						
						my $n_new_rule = $self->fire_rule($$facts_ref, $rule, $constraints);
						$self->debug_print( ($n_new_rule ? 33 : 34) , "   ".join("\t", $n_lhs, rule_to_string($rule)) ) ;
						$n_new += $n_new_rule ;
						$visited_rules{$match_id} = 1;
					}
				} else {
# 					$self->debug_print( 30, "   ".join("\t", join("/", scalar( @{ $$constraints{'matched'} } ), scalar( @{ $$constraints{'unmatched'} } )), rule_to_string($rule)) ) 
# 						if scalar( @{ $$constraints{'unmatched'} } ) == 1;
				}
			}
			last if $n_new > 0;
		}
		
		TO_NEXT_CYCLE:
		$self->debug_print(36, $$facts_ref->get_fact_tag_strings() ) ;
	} while ($n_new != 0);
	
# 	$facts->clean_tags();
	
	return $$facts_ref; 	# $facts->get_fact_tag_hash();
}

our %__Rule_load_proper_cache ;
our %__Rule_load_proper_cache_long ;

sub __Rule_load_proper($$$$) {
	my ($self, $order, $prio, $rule_str) = @_;
	
	return undef if exists $self->{'rules_str'}{$rule_str};
	
	my ($lhs, $rhs) = split /\s*=>\s*/, $rule_str;
	
	return undef if ! defined $rhs;

	my @lhs ;
	my @lhs_neg ;
	my @lhs_assert ; # false assertion for conditional constraint satisfication. Syntax "*LHS1"
	
	if ( exists $__Rule_load_proper_cache_long{$lhs} ) { # Two level caching to improve performance.
		my $c = Storable::thaw( $__Rule_load_proper_cache_long{$lhs} ) ;
		@lhs_assert = @{ $$c{a} };
		@lhs_neg = @{ $$c{n} };
		@lhs = @{ $$c{e} };
	} else {
		for my $lhse ( split /\s*;\s*/, $lhs ) {
			my $lhse0 = $lhse ;
			if ( exists $__Rule_load_proper_cache{$lhse0} ) {
				my $c = Storable::thaw( $__Rule_load_proper_cache{$lhse0} ) ;
				push @lhs_assert, $$c{a};
				push @lhs_neg, $$c{n};
				push @lhs, $$c{e};
			} else {
				$lhse =~ s/^\s*//; 
				$lhse =~ s/\s*$//; 
				push @lhs_assert, ( my $lhsa = ( ( $lhse =~ s/^\s*\*\s*// ) ? 1 : 0) ); 
				push @lhs_neg, ( my $lhsn = ( ( $lhse =~ s/^\s*NOT\b\s*// ) ? 1 : 0) ) ; 
				push @lhs, $lhse;
				$__Rule_load_proper_cache{$lhse0} = Storable::freeze( { a => $lhsa, n => $lhsn, e => $lhse } ) ;
			}
		}
		$__Rule_load_proper_cache_long{$lhs} = Storable::freeze( { a => \@lhs_assert, n => \@lhs_neg, e => \@lhs } ) ;
	}

	my @rhs ;
	my @rhs_neg ;
	
	if ( exists $__Rule_load_proper_cache_long{$rhs} ) { # Two level caching to improve performance.
		my $c = Storable::thaw( $__Rule_load_proper_cache_long{$rhs} ) ;
		@rhs_neg = @{ $$c{n} };
		@rhs = @{ $$c{e} };
	} else {
		my $rhs0 = $rhs;
		while ( length($rhs) ) {
			$rhs =~ s/^\s*//;
			if ( $rhs =~ /^(.+?(?:\[.+?\])?)\s*(?:;|$)/ ) {
				$rhs = $';
				
				my $rhse0 = my $rhse = $1;
				if ( exists $__Rule_load_proper_cache{$rhse0} ) {
					my $c = Storable::thaw( $__Rule_load_proper_cache{$rhse0} ) ;
					push @rhs_neg, $$c{n};
					push @rhs, $$c{e};
				} else {
					push @rhs_neg, ( my $rhsn = ( ( $rhse =~ s/^\s*NOT\b\s*// ) ? 1 : 0) ); 
					push @rhs, $rhse;
					$__Rule_load_proper_cache{$rhse0} = Storable::freeze( { a => 0, n => $rhsn, e => $rhse } ) ;
				}
			}
		}
		$__Rule_load_proper_cache_long{$rhs0} = Storable::freeze( { a => [ 0 x scalar(@rhs) ], n => \@rhs_neg, e => \@rhs } ) ;
	}
	
	my %a = ('order' => $order, 'prio' => $prio, 'lhs' => \@lhs, 'lhs_neg' => \@lhs_neg, 'lhs_assert' => \@lhs_assert, 'rhs' => \@rhs, 'rhs_neg' => \@rhs_neg);
	
	$self->{'rules_str'}{$rule_str}++;
	
	return \%a;
}


sub load { # Load if-then rules from strings
	my ($self, @rules_strs) = @_;
	my @rules_spec;
	
	my %disjunctive_lhs;

	my $order = 1;

	for my $rule_str (@rules_strs) {
		$rule_str  =~ s/\s*#.*//; # remove comments
		my $prio = 0; # defines rule priority;
		my $prefix = ( $rule_str =~ s/^\s*\[(.+?)\]\s*// ) ? $1 : undef;
		
		if ( defined $prefix ) {
			my %attrib = map { my ($a, $b) = split /\s*:\s*/, $_; $a => $b } split /\s*;\s*/, $prefix ;
			$prio = $attrib{'prio'} if exists $attrib{'prio'} ;
		}
		
		my $a = $self->__Rule_load_proper($order, $prio, $rule_str);
		next if ! defined $a;
		++$order ;
		push @rules_spec, $a ; # compiled rules
		
		$disjunctive_lhs{$_} = $a for grep { /^\(.*?\s+OR\s+.*\)$/ } @{ $a->{lhs} } ; # Now can handle simple disjunctions (one level only) FIXME
	}

	for my $disj_lhs (sort { $disjunctive_lhs{$a}{'order'} <=> $disjunctive_lhs{$b}{'order'} } keys %disjunctive_lhs) {
		my ($disj_rhs, $disj_rhs_orig) = ($disj_lhs, $disj_lhs);
		$disj_rhs =~ s/^\((.+)\)$/$1/;
		my @disj_lhs = split /\s+OR\s+/, $disj_rhs ;
		
		my $prio = ( $disjunctive_lhs{$disj_lhs}{'prio'} - 1 ); # Handling of rules should take precedence before the actual rules.
		
		for my $disj_lhs (@disj_lhs) { 
			my $rule_str = "$disj_lhs => $disj_rhs_orig";
			
			my $a = $self->__Rule_load_proper($order, $prio, $rule_str);
			next if ! defined $a;
			
			++$order ;
			push @rules_spec, $a; # compiled rules
		}
	}
	
	my $n_rules = scalar @rules_spec;
	$self->debug_print(37, "Added $n_rules rules.");
	push @{ $self->{'rules'} }, @rules_spec;
	$self->{'f_needs_reindex'} = 1;
}


sub define_dyn_rule {
	my ($self, $rule_name, $code_ref) = @_;
	my $n_drules = scalar keys %{ $self->{'dyn_rules'} };
	my $drule_id = "dr-". $n_drules ;
	
	$self->{'dyn_rules'}{$drule_id}{'name'} = $rule_name;
	$self->{'dyn_rules'}{$drule_id}{'code'} = $code_ref;
}

sub facts_to_graphviz {
	my @facts = @_;
	our $last_fact_id = 0;
	our %fact_to_id;
	our %id_to_fact;
	
	sub get_fact_id {
		my $f = shift;
		if ( ! exists $fact_to_id{$f} ) {
			my $id = sprintf("F%04d", $last_fact_id);
			$fact_to_id{$f} = $id;
			$id_to_fact{$id} = $f;
			++$last_fact_id;
		}
		return $fact_to_id{$f} ;
	}
	
	for my $line (@facts) {
		chomp $line;
		my $fact = $line;
		my @tags ;
		if ( $fact =~ /(.*)\s*\[(.+)\]$/ ) {
			$fact = $1;
			@tags = split /\s*;\s*/, $2;
		}
		print join(" ", map { "[$_]" } (get_fact_id ($fact), @tags) );
		print "\n";
	}
}

##################################################################################
# OLD TEST CODE:
# my $ruleset = Rules->new( {debug_level=>1});
# $ruleset->load('A => B; C [tag1]');
# $ruleset->load('A; NOT B => E');
# $ruleset->load('A; B => D [tag2]');
# $ruleset->load('A; (C OR E) => D [tag3]');
# $ruleset->load('A; B => D [tag2]');
# $ruleset->load('D => A [tag3]');
# my $facts = Facts->new('A', 'E') ;
# my %output = $ruleset->run( \$facts ) ->get_fact_tag_hash();
# # my %output = $ruleset->run( \( Facts->new('A', 'E') ) ) ->get_fact_tag_hash();
# $ruleset->print;
# print map { "$_\t[$output{$_}]\n" } sort keys %output;
# 
# 1;
}


{
package Rules::Modules ;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 0.8 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}



sub new {
    my ($class, $params) = @_;
    my $self = {};
    
    $self->{'_params'} = $params ;
    $self->{'modules'} = {};
    
    bless($self, $class);
    
    return $self;
}


sub add_module {
	my ($self, $ruleset_name) = @_;
# 	my $cnt = $self->{'modules'}
	$self->{'modules'}{$ruleset_name} = Rules->new( \%{ $self->{'_params'} } )
		if ! exists $self->{'modules'}{$ruleset_name};
	$self->{'_latest_module'} = $ruleset_name;
	return $self->{'modules'}{$ruleset_name} ;
}

sub load_rules {
	my ($self, $ruleset_name, @rules) = @_;
	$ruleset_name = $self->{'_latest_module'} if ! defined $ruleset_name or ! length $ruleset_name;
	my $ruleset = $self->add_module($ruleset_name);
	
	$self->debug_print(37, "Loading ".scalar(@rules)." into [$ruleset_name]: ");
	$ruleset->load( @rules );
}

sub print {
	my $self = shift;
	for my $ruleset_name (sort keys %{ $self->{'modules'} }) {
		$self->{'modules'}{$ruleset_name}->print();
	}
}

sub run {
	my ($self, @facts) = @_;

	$self->{'cache'} = (); # clear the runtime cache, used for dynamic storage of runtime private facts
	$self->{'debug_output'} = '';
	
	my $debug_output ;
	
	my $facts = Facts->new(@facts) ;
	
	for my $ruleset_name (sort keys %{ $self->{'modules'} }) {

		$self->{'modules'}{$ruleset_name}->debug_print(37, "Running module: $ruleset_name");
		
		$self->{'modules'}{$ruleset_name}->run(\$facts);

		$debug_output .= $self->{'modules'}{$ruleset_name}{'debug_output'};
	}

	$self->{'debug_output'} = $debug_output ;
	
	return $facts->get_fact_tag_strings();
}


1;
}

# my $rule = Rules->new();
# 
# my @fact_patts = ('is-a:$X:?Y; has:?X => has:?Y');
# 
# $rs->define_dyn_rule( 'process_fact_patterns',  sub { $_=shift; 
# 		for my $patt (@fact_patts) {
# 			return () if ! /^is-a:(\S+?):(\S+?)/ and ! /(\S+?)/ ;
# 		}

# );

# my $r = Rules->new(); 
# $r->load("[prio:0] A=>B", "[prio:-1] B=>C", "[prio:0] C=>D"); 
# $r->print; 
# my %l = $r->run("A"); 
# print map {"$_\t$l{$_}\n"} sort keys %l;
# print $r->{'debug_output'} ;

# my $r = Rules->new( { run_rules_only_once => 1}  ); 
# my $f = Facts->new(); 
# $r->load("0 => A", "A => NOT A; B", "B => NOT B; C", "C => NOT C; D", "D => NOT D; A"  );  
# $r->print; 
# $f->assert("A");
# my %l = $r->run(\$f)->get_fact_tag_hash(); 
# print map {"$_\t$l{$_}\n"} sort keys %l;
# print $r->{'debug_output'} ;

1;
