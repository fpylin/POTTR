#!/usr/bin/perl
##############################################################################
#
# TSV.pm - Perl class for handling Tab-separated value dataframes
# 
# Copyright 2015-2019, Frank Lin
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

package TSV;
# use Devel::StackTrace;
# my $trace = Devel::StackTrace->new;

use strict;
use warnings;

BEGIN {
	use Exporter   ();
	our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = 1.00;
	$VERSION     = sprintf "%d.%03d", q$Revision: 1.1 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
	@EXPORT      = qw();
	%EXPORT_TAGS = ();
	@EXPORT_OK   = qw();
    }

###############################################################################
sub v_safe { my $x = shift; return ( defined $x ? $x : "" ); }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub shuffle(\@) { # Fisher–Yates shuffle Algorithm
	my $l = shift;
	my $n = scalar(@$l) ;
	for(my $i =$n-1; $i >=1; --$i ) {
		my $j = int( rand( $i ) );
		my $s = ${$l}[$i];
		${$l}[$i] = ${$l}[$j];
		${$l}[$j] = $s;
	}
	return @$l;
}
sub max { return undef if (! scalar(@_) );
	my $v = shift;
	for (@_) { $v = $_ if ($_ > $v) ; }
	return $v;
}
sub argmax(\%) { my %a = %{ $_[0] }; my $x = undef;
	for (keys %a) { $x = $_ if ((! defined $x) || ($a{$_} > $a{$x})); }
	return $x;
}


###############################################################################
sub new {
    my ($class, $filename, $field_ref) = @_;
    
    my $self = { 'filename' => $filename };
    
    bless($self, $class);
    
    if ( defined($filename) and length($filename) ) {
        $self->load($filename, $field_ref) ;
    }
    
    return $self;
}

###############################################################################
sub new_from_data {
    my ($class, @lines) = @_;
    
    my $self = { 'filename' => '' };
    
    bless($self, $class);
    
    $self->import_data( @lines );
    
    return $self;
}


###############################################################################
sub load {
    my $self = shift;
    my $filename = shift;
    my $field_ref = shift;
    
#    print STDERR "$filename: ";
    
    my @lines = file($filename);
    
    chomp for @lines ;
    
    unshift @lines, join("\t", @{ $field_ref }) if defined $field_ref ;
    
    return $self->import_data( @lines );
}


###############################################################################
sub import_data {
    my $self = shift;    
    my @lines = @_;
    
    chomp for @lines ;

    my @clean;
    $self->{'data'} = \@clean;
    
    my $header = shift @lines ;
    my @body = @lines ;
    
#    print STDERR join(", ", caller()) if ! defined $header;
    
    my @fields = split /\t/, $header;
    my %fields ;
    for my $f (@fields) {
        if ( ! exists $fields{$f} ) {
            $fields{$f} ++;
            next;
        }
        
        for (my $cnt=2; $cnt< 9999; ++$cnt) {
            my $varname = "$f.$cnt";
            next if exists $fields{$varname} ;
            $fields{$varname} ++;
            $f = $varname;
            last;
        }
    }
    
    $self->{'fields'} = \@fields;
    
    for my $line (@body) {
        my %a ;
        my @values = split /\t/, $line ;
        for my $i ( 0 .. $#fields ) {
#             next if ! defined $values[$i] ;
            $a{ $fields[$i] } = $values[$i] ;
        }
        push @{ $self->{'data'} }, \%a;
    }
    
#    print STDERR $self->n_fields()." col x ".$self->n_rows()." rows imported.\n";
    
    return scalar(@body);
}

###############################################################################
sub set_fields {
	my $self = shift;
    my @fields = @_;
    $self->{'fields'} = \@fields ;
}

###############################################################################
sub push_rows {
	my $self = shift;
	my @rowrefs = @_;
	for my $rr (@rowrefs) {
		die if ref($rr) ne 'HASH';
		push @{ $self->{'data'} }, $rr;
	}
}

###############################################################################
sub get_column {
	my $self = shift;
	my $field = shift;
	my @values = map { ${$_}{$field} } @{ $self->{'data'} };
# 	pop @values while scalar @values and ! defined $values[$#values];
	return @values ;
}

###############################################################################
sub import_column {
	my $self  = shift;
	my $field = shift;
	my @data  = @{ $_[0] } ;
	
	push @{ $self->{'fields'} }, $field;
	my $N = undef;
	
 	print STDERR "TSV::import_column - ".scalar(@data)."\n";
# 	print STDERR join("\n", map { "$_\t$data[$_]" } ( 0..$#data ) )."\n";
	
	if ( ! defined $self->{'data'} ) {
		$self->{'data'} = ();
		$N = scalar(@data);
	}
	else {
		$N = $self->n_rows() ;
		die "Number of rows ($N) mismatch: imported data has ".scalar(@data)." rows.\n" if $self->n_rows() != scalar(@data);
	}
	
	for ( my $i=0; $i<$N; ++$i )  {
		$self->{'data'}[$i]{$field} = $data[$i];
	}
}

###############################################################################
sub save_as {
    my $self = shift;
    my $filename = shift;
    my $params = shift;
    die "No filename specified" if ! defined $filename;
    
    my @fields = @{ $self->{'fields'} };

    open FOUT, ">$filename" or die "Unable to open $filename for writing.\n";
    print FOUT join("\t", @fields)."\n" if ( (! exists $self->{'do-not-save-header'}) and (! defined $params or ! exists $$params{'no-header'}) );
    for my $row ( @{ $self->{'data'} } ) {
        my %a = %{ $row };
        print FOUT join("\t", map { v_safe($a{$_}) } @fields)."\n";
    }
    close FOUT;
    
    print STDERR "Writing ".$self->n_fields()." col x ".$self->n_rows()." rows to file $filename\n";
}

###############################################################################
sub save_R_table {
    my $self = shift;
    my $filename = shift;
    die "No filename specified" if ! defined $filename;
    
    my @fields = @{ $self->{'fields'} };

    open FOUT, ">$filename" or die "Unable to open $filename for writing.\n";
    print FOUT join("\t", @fields)."\n";
    my $cnt = 1;
    for my $row ( @{ $self->{'data'} } ) {
        my %a = %{ $row };
        print FOUT join("\t", $cnt, map { $a{$_} } @fields)."\n";
        ++$cnt;
    }
    close FOUT;
    
    print STDERR "Writing ".$self->n_fields()." col x ".$self->n_rows()." rows to file $filename.\n";
}

###############################################################################
sub to_string {
    my $self = shift;
    my @output;
    
    my @fields = @{ $self->{'fields'} };
    
    push @output, join("\t", @fields)."\n";
#     print STDERR "TSV->get(): ".join(", ", @fields)."\n";
    
    for my $row ( @{ $self->{'data'} } ) {
        my %a = %{ $row };
        my @z ; warn "TSV->get(): Undefined field: ".join(", ", @z)."\n" if ( @z = grep { ! defined $a{$_} } @fields );
        push @output,  join("\t", map { defined $a{$_} ? $a{$_} : 'NA' } @fields)."\n";
    }
    
#     print STDERR "TSV->get(): output\n".join("", @output)."\n";
    return @output;
}


###############################################################################
sub ith_field {
    my $self = shift;
    my $i = shift;
    return @{ $self->{'fields'} }[$i];
}

###############################################################################
sub n_fields {
    my $self = shift;
    return scalar( @{ $self->{'fields'} } );
}

###############################################################################
sub last_field {
    my $self = shift;
    my $n = scalar( @{ $self->{'fields'} } );
    return ${ $self->{'fields'} }[$n-1];
}

###############################################################################
sub guess_class_label {
    my $self = shift;
    return 'class' if grep { $_ eq 'class' } @{ $self->{'fields'} } ;
    return $self->last_field;
}

###############################################################################
sub n_rows {
    my $self = shift;
    return scalar( @{ $self->{'data'} } );
}


###############################################################################
sub levels {
    my $self = shift;
    my $field = shift;
    my @x = map { ${$_}{$field} } @{ $self->{'data'} };
    my @u = uniq(@x);
    return @u;
}

###############################################################################
sub foreach_row {
	my $self = shift;
	my $sub = shift;
	my @args = @_;
	for my $row ( @{ $self->{'data'} } ) {
		$sub->($row, @args);
	}
}

###############################################################################
sub enumerate_feature_levels {
    my $self = shift;
    my $field = shift;
    my $feature_level_str = '';
    my $regex_T = '(?i:T|TRUE|P|Present|Yes)';
    my $regex_F = '(?i:F|FALSE|F|Absent|No|MISSING)';
	my %histogram = $self->histogram($field);
	my @levels = keys %histogram ;
	my @non_na_levels = grep { ! /^(?:NA|\?)$/ } @levels ;
	my @non_numeric_levels = grep { ! /^(?:[\+\-]?[0-9]+(?:\.?[0-9]+)?(?:[Ee][\+\-]?[0-9]+)?)$/ } @non_na_levels ;
	
	if ( scalar(@non_numeric_levels) == 0 ) { # is numeric
		$feature_level_str = 'NUMERIC';
	} else {
		my $f_has_na = ( scalar(@levels) != scalar(@non_na_levels) ) ? 1 : 0;
		my $base = argmax(%histogram); # Choose mode as the base
		my %weights = %histogram ;
		if ( scalar( grep {/^(?:regex_F)$/ } @non_na_levels ) + scalar( grep {/^(?:regex_T)$/ } @non_na_levels ) == 2 ) {
			my ($base) = ( grep {/^regex_F$/ } @non_na_levels);
			$weights{$base} = max(values %histogram) + 1;
        }
        my @order = sort { $weights{$b} <=> $weights{$a} } @non_na_levels;
        $feature_level_str = join("\t", @order);
    }
    $self->{'feature_levels'}{$field} = $feature_level_str;
}

sub enumerate_all_feature_levels {
    my $self = shift;
    for my $f ( @{ $self->{'fields'} } ) {
		$self->enumerate_feature_levels($f);
    }
}

###############################################################################
sub histogram {
    my $self = shift;
    my $field = shift;
    my %cnt ;
    for my $row ( @{ $self->{'data'} } ) {
        my $v = ${$row}{$field};
        if ( ! defined $v ) {
# 			warn "Value for field ''$field'' not defined" 
		} else {
			$cnt{ $v }++;
		}
    }
    return %cnt;
}

###############################################################################
sub dummify { # dummify everything except the arguments
    my $self = shift;
    my %fields_except = map { $_=>1 } @_;
    
    my @fields = @{ $self->{'fields'} };
    my @fields_new;
    my @fields_to_delete;
    
    my $regex_T = '(?i:T|TRUE|P|Present|Yes)';
    my $regex_F = '(?i:F|FALSE|F|Absent|No)';
    for my $f ( @fields ) {
		my %histogram = $self->histogram($f);
		my @levels = keys %histogram ;
		my @non_na_levels = grep { ! /^(?:NA|\?)$/ } @levels ;
		my @non_numeric_levels = grep { ! /^(?:[\+\-]?[0-9]+(?:\.?[0-9]+)?(?:[Ee][\+\-]?[0-9]+)?)$/ } @non_na_levels ;		
		
		my $f_has_na = ( scalar(@levels) != scalar(@non_na_levels) ) ? 1 : 0;
		
#         print "$f:", join(" ", map { "[$_]"} sort @non_numeric_levels)."\n";
		
        if ( ( exists $fields_except{$f} ) or ( scalar(@non_numeric_levels) == 0 ) ) {
			push @fields_new, $f;
			next;
		}
# 		print STDERR "Dummifying field ''$f''\n";
		
		my $base =  argmax(%histogram); # Choose mode as the base
		if ( scalar( grep {/^(?:regex_F)$/ } @non_numeric_levels ) + scalar( grep {/^(?:regex_T)$/ } @non_numeric_levels ) == 2 ) {
			($base) = ( grep {/^regex_F$/ } @non_numeric_levels);
		}
		push @fields_to_delete, $f;
		for my $l (@non_numeric_levels) {
			next if $l eq $base;
			my $varname = "$f.$l";
			push @fields_new, $varname;
			for my $row ( @{ $self->{'data'} } ) {
				${$row}{$varname} = ( ${$row}{$f} eq $l ) ? 1 : 0;
			}
		}
	}
	$self->{'fields'} = \@fields_new;
}

###############################################################################
sub fields {
    my $self = shift;    
    return @{ $self->{'fields'} } ;        
}

###############################################################################
sub remove {
    my $self = shift;    
    my @fields_to_delete = @_;
    my $fields_to_delete_regex = join('|', @fields_to_delete);
    my @fields_new = grep { ! /^(?:$fields_to_delete_regex)$/ } @{ $self->{'fields'} } ;
    $self->{'fields'} = \@fields_new ;
    my %fields_new = map { $_=>1 } @fields_new ;
    for my $f ( keys %{ $self->{'feature_levels'} } ) {
		delete $self->{'feature_levels'}{$f} if ! exists $fields_new{$f};
    }
    return scalar(@{ $self->{'fields'} });
}

###############################################################################
sub remove_rows_by_criteria {
    my $self = shift;
    my $code = $_[0];
    my @new_data;
    for my $row ( @{ $self->{'data'} } ) {
        my %a = %{ $row } ;
        next if $code->(\%a);
        push @new_data, $row;
    }
    $self->{'data'} = \@new_data ;
}

###############################################################################
sub translate {
    my $self = shift;
    my $field = shift;
    my $xlate_tab_ref = shift;
    
    for my $row ( @{ $self->{'data'} } ) {
		my $v = ${$row}{$field};
		${$row}{$field} = ${$xlate_tab_ref}{$v} if exists ${$xlate_tab_ref}{$v};
	}
}


###############################################################################
sub is_numeric {
    my $self = shift;
    my $field = shift;
    my @levels = grep { ! /^NA$|^\?$|^\s*$/ } $self->levels($field);
    my @n_match = grep { /^[+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?$/ } @levels;
    return undef if scalar(@levels) == 0;
    return (scalar(@n_match) == scalar(@levels));
}

###############################################################################
sub select { # Select by fields
    my $self = shift;
    my @wanted_fields = @_;
    
    my $TSV_new = new TSV;
    
    $TSV_new->{'fields'} = \@wanted_fields;

    my @data;
    
    for my $row ( @{ $self->{'data'} } ) {
        my %a = map { $_ => ${$row}{$_} } @wanted_fields;
        push @data, \%a;
    }
    $TSV_new->{'data'} = \@data;
    
    return $TSV_new;
}

###############################################################################
sub subsample { # Select by row numbers
    my $self = shift;
    my @irows = @_;

    my $TSV_new = new TSV;
    
    $TSV_new->{'fields'} = $self->{'fields'} ;
    
    my @a = @{ $self->{'data'} }[ @irows ] ;
    $TSV_new->{'data'} = \@a;
    
    return $TSV_new ;
    }

    
###############################################################################
sub grep_rows{ 
    my $self = shift;
    my $code_ref = $_[0];

    my $TSV_new = new TSV;
    
    $TSV_new->{'fields'} = $self->{'fields'} ;
    
    my @a = grep { $code_ref->($_) } @{ $self->{'data'} };
    
    $TSV_new->{'data'} = \@a;
    
    return $TSV_new ;
    }

###############################################################################
sub shuffle_rows { # shuffling data
    my $self = shift;
    
    @{ $self->{'data'} } = shuffle( @{ $self->{'data'} } );
    }

###############################################################################
sub stratify_by { # stratify data by field
	my $self = shift;
	my $f_strata = shift;
	my @levels = $self->levels($f_strata);

	my $N = scalar( @{ $self->{'data'} });
	my %by_class ;
	
	for my $row ( @{ $self->{'data'} } ){
		my $class_label = ${$row}{$f_strata};
		push @{ $by_class{$class_label} }, $row ;
	}
	
	sub sum { my $sum = 0; $sum += $_ for @_; return $sum; }

	my %expected ;
	for my $cl (@levels) {
		$expected{$cl} = ( scalar(@{ $by_class{$cl} }) / $N );
	}
	
	@levels = sort { scalar(@{ $by_class{$b} }) <=> scalar(@{ $by_class{$a} }) } @levels ;

	my @new_data;
	
	while ( (my $s = sum( map { scalar( @{ $by_class{$_} } ) } @levels )) > 0 ) {
		for my $cl (@levels) {
			my $Oclass= scalar(@{$by_class{$cl}}) / $s;
			my $Eclass= $expected{$cl};
			if ( $Oclass >= $Eclass ) {
				my $o = shift @{ $by_class{$cl} };
				push @new_data, $o;
				last;
			}
		}
	}
	$self->{'data'} = \@new_data;
}
    
    
###############################################################################
sub rename { # rename field
    my $self = shift;
    my $x = shift;
    my $y = shift;
    
    my $nfields = scalar @{ $self->{'fields'} };
    for ( my $i=0; $i<$nfields; ++$i ) {
        next if ( $self->{'fields'}[$i] ne $x ) ;
        $self->{'fields'}[$i] = $y;
        for ( @{ $self->{'data'} }) {
            ${$_}{$y} = ${$_}{$x};
            delete ${$_}{$x};
        }
    }
}


###############################################################################
sub try_quote_WEKA {
    my $x = shift;
    return '?' if ! defined $x or $x =~ /^\s*$/ or $x eq 'NA';
    return "\"$x\"" if $x =~ /[\s\+]/;
    return $x;
}


###############################################################################
sub space_safe {
	my $x = shift;
	$x =~ s/ /_/g;
	return $x ;
}

###############################################################################
sub to_arff {
    my $self = shift;
	my @output;
	my @fields = @{ $self->{'fields'} };
	
	push @output, "\@relation relation1\n";
	
	$self->enumerate_feature_levels() if ! defined $self->{'feature_levels'};

	for my $f (@fields) {
		my @levels;
        if ( exists $self->{'feature_levels'}{$f} ) {
			@levels = split /\t/, $self->{'feature_levels'}{$f} ;
# 			print STDERR "\e[1;32m$f".join("!", @levels)."\e[0m\n";
		} else {
			@levels = $self->levels($f);
# 			print STDERR "\e[1;34m$f".join("!", @levels)."\e[0m\n";
		}
		
		if ( ( scalar(@levels) && ($levels[0] eq 'NUMERIC') ) || $self->is_numeric($f) ) {
            push @output, "\@attribute ".space_safe($f)." NUMERIC\n"; 
            next;
        }
        
		# Categorical data
# 		print STDERR "\e[1;36m".join("!", @levels)."\e[0m\n";
		push @output, "\@attribute ".space_safe($f)." {".join(", ", grep { $_ !~ /^(?:\?|NA)$/ } map { try_quote_WEKA($_) } @levels)."}\n"; 
	}
        
	push @output, "\@data\n";
	for my $row ( @{ $self->{'data'} } ) {
        my %a = %{ $row };
        push @output, join(",", map { try_quote_WEKA($a{$_}) } @fields)."\n";
	}
	
	return @output;
}

###############################################################################
sub index_by{
    my $self = shift;
    my $field = shift;
    my %a ;
    for my $row ( @{ $self->{'data'} } ) {
		warn "Field $$row{$field} already exists\n" if exists $a{ $$row{$field} };
		$a{ $$row{$field} } = $row;
    }
    return %a;
}


END { }       # module clean-up code here (global destructor)

1;  # don't forget to return a true value from the file
