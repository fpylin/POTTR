#!/usr/bin/perl
##############################################################################
#
# Config.pm - precision oncology trial and therapy recommender
# 
# Configuration module
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

package POTTRConfig;

use strict;
use warnings;

use POSIX;
use Cwd;

use lib '.';

sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 1.1 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}

our %data;
our @predefined_rules;
our $base_rel_path ;

my $f_initiailised ;

sub import {
	my $self = shift;
	$base_rel_path = shift if scalar @_;
# 	print STDERR "base_path = $base_rel_path. called from ".(join ":", caller)."\n" if defined $base_rel_path;
# 	print STDERR "No base_path defined. called from ".(join ":", caller)."\n"  if ! defined $base_rel_path;
}

sub load {
	my $config_file = shift;
# 	print STDERR "Config file $config_file loaded.\n";
	for ( file($config_file) ) {
		s/#.*|//g;
		s/^\s*|\s*$//g;
		if ( /^set +(.+)/i ) {
			my $x = $1;
			if ( $x =~ /=>/ ) { # is a rule
				push @predefined_rules, $x;
			} else {
				push @predefined_rules, '(initial-fact) => '.$x;
			}
			
		}
		
		if ( my ($lhs, $rhs) = split /\s*=\s*/, $_ ) {
			push @{ $data{$lhs} }, $rhs;
		}
	}
	$f_initiailised = 1;
}

sub get {
	&ON_DEMAND_INIT;
	my $x = shift;
	return () if ! exists $data{$x} ;
	return @{ $data{$x} };
}

sub get_paths {
	&ON_DEMAND_INIT;
	my $type = shift;
	my $x = shift;

	die "FATAL: $type directory not defined.\n".(join(" ", keys %data))."\n" if ! exists $data{$type.'_dir'};
	my @dirs = $data{$type.'_dir'};
	my $dir  = $data{$type.'_dir'}[0];
	
	$dir = "$base_rel_path/$dir" if ( $dir !~ /^\// ) and ( defined $base_rel_path ) and ( -d "$base_rel_path/$dir" ) ;
	
	if ( ! -d $dir ) {
		die "FATAL: [$type] directory ``$dir'' does not exist at ".join(":", caller)."\n" ;
	}
	
	return () if ! exists $data{$x} ;
	return map { "$dir/$_" } @{ $data{$x} };
}

sub get_first_path {
	&ON_DEMAND_INIT;
	my $type = shift;
	my $x = shift;
	my @retval = get_paths($type, $x);
	return $retval[0];
}

sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;
	
	load( "pottr_config" );
}

1;
