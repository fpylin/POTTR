#!/usr/bin/perl
##############################################################################
#
# DOID.pm - Disease Ontology interface
#
# This file is part of Oncology Treatment and Trials Recommender (OTTR) and 
# Precision OTTR (POTTR).
#
# Copyright 2019, Frank Lin & Kinghorn Centre for Clinical Genomics, 
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
#
# Interm use of Disease Ontology -- to upgrade to a precision oncology-specific 
# ontology in future version
#
#
#

use strict;
use warnings;

package DOID;

use lib '.';
use POTTRConfig;

sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }
sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub mtime { my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($_[0]); return $mtime; }

##################################################################################

our $f_initiailised = undef;

our @EXPORT;
our @EXPORT_OK;

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = sprintf "%d.%03d", q$Revision: 1.1 $ =~ /(\d+)/g;
    @ISA         = qw(Exporter);
    @EXPORT      = qw( %DO_name %DO_id %DO_is_a $DO_regex &DO_match_catype );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = qw();
}

our %DO_name;
our %DO_is_a;
our %DO_id;

my @POTTR_parts = (q|[Term]
id: POTTR:0001
name: solid tumour
def: "solid tumours"
subset: DO_cancer_slim
synonym: "Solid tumour" EXACT []
synonym: "tumors" EXACT []
synonym: "Tumors" EXACT []
synonym: "Tumor" EXACT []
synonym: "advanced cancer" EXACT []
synonym: "advanced solid tumour" EXACT []
synonym: "Advanced Tumors" EXACT []
synonym: "pan cancer" EXACT []
synonym: "other tumor types" EXACT []
synonym: "Neoplasm" EXACT []
synonym: "Neoplasms" EXACT []
is_a: DOID:162 ! Neoplasm|, 
q|[Term]
id: POTTR:0002
name: cancer of unknown primary
def: "cancer of unknown primary"
subset: DO_cancer_slim
synonym: "carcinoma of unknown primary" EXACT []
is_a: DOID:162 ! Neoplasm|, 
q|[Term]
id: POTTR:1000
name: Diffuse Intrinsic Pontine Glioma
def: "Diffuse Intrinsic Pontine Glioma"
subset: DO_cancer_slim
synonym: "DIPG" EXACT []
is_a: DOID:162 ! Neoplasm|, 
q|[Term]
id: POTTR:1001
name: triple-negative breast cancer
def: "Triple-negative breast cancer"
subset: DO_cancer_slim
synonym: "TNBC" EXACT []
is_a: DOID:1612 ! Breast cancer|, 
q|[Term]
id: POTTR:1002
name: GastroEsophageal Cancer
def: "GastroEsophageal Cancer"
subset: DO_cancer_slim
synonym: "Gastrooesophageal Cancer" EXACT[]
synonym: "Esophagogastric Cancer" EXACT[]
synonym: "Esophagogastric Cancers" EXACT[]
synonym: "Gastroesophageal junction adenocarcinoma" EXACT[]
synonym: "Gastroesophageal junction cancer" EXACT[]
synonym: "Gastroesophageal junction carcinoma" EXACT[]|, 
q|[Term]
id: POTTR:1003
name: Chronic myelogenous leukaemia
def: "Chronic Myelogenous Leukaemia"
subset: DO_cancer_slim
synonym: "Chronic Myelogenous Leukaemia" EXACT[]
synonym: "Chronic Myelogenous Leukemia" EXACT[]
synonym: "CML" EXACT[]
synonym: "chronic myeloid leukaemia" EXACT[]
is_a: POTTR:1007 ! Liquid cancer|, 
q|[Term]
id: POTTR:1004
name: Acute lymphoblastic leukaemia
def: "Acute Lymphoblastic Leukaemia"
subset: DO_cancer_slim
synonym: "Acute lymphoblastic leukemia" EXACT[]|, 
q|[Term]
id: POTTR:1005
name: Colon cancer
def: "Colon Cancer"
subset: DO_cancer_slim
is_a: DOID:0080199 ! Colorectal carcinoma|,
q|[Term]
id: POTTR:1006
name: Rectal cancer
def: "Rectal Cancer"
subset: DO_cancer_slim
is_a: DOID:0080199 ! Colorectal carcinoma|,
q|[Term]
id: POTTR:1007
name: Liquid cancer
def: "Liquid cancer"
synonym: "Liquid cancers"
subset: DO_cancer_slim
is_a: DOID:162 ! Neoplasm|, 
q|[Term]
id: POTTR:1008
name: Glioma
def: "Glioma"
synonym: "Gliomas"
subset: DO_cancer_slim
is_a: POTTR:0001 ! Solid tumour|,
q|[Term]
id: POTTR:1009
name: Inflammatory myofibroblastic tumour
def: "Inflammatory myofibroblastic tumour"
synonym: "Inflammatory myofibroblastic tumor"
synonym: "Inflammatory myofibroblastic tumour"
synonym: "Inflammatory myofibroblastic tumors"
synonym: "Inflammatory myofibroblastic tumours"
subset: DO_cancer_slim
is_a: POTTR:0001 ! Solid tumour|,
q|[Term]
id: POTTR:1010
name: Aggressive fibromatosis
def: "Aggressive fibromatosis"
synonym: "Desmoid"
synonym: "Desmoid, NOS"
subset: DO_cancer_slim
is_a: POTTR:0001 ! Solid tumour|,
);

our @synonym_cancer = ('cancer', 'tumor', 'tumour', 'carcinoma', 'adenocarcinoma', 'squamous cell carcinoma', 'neoplasm', 'malignancy', 'malignancies');
our %synonym_map = (
	'-' => ' ',
	'-' => '',
	'b-' => 'b-cell ',
    'head and neck squamous cell carcinoma' => 'Squamous Cell Carcinoma of the Head and Neck|SCCHN|HNSCC|Locally advanced Squamous Cell Carcinoma of the Head and Neck|Neoplasms, Head and Neck',
    'hepatocellular carcinoma' => 'Carcinoma, Hepatocellular',
# 	'gastroesophageal' => 'esophagogastric',
#     'gastro-esophageal' => 'gastroesophageal|gastro-oesophageal junction|Gastro-Oesophageal',    
# 	'esophagus' => 'oesophagus|oesophageal|esophageal|Oesophageal',
	'skin' => 'cutaneous',
	'cutaneous' => 'skin',
	'squamous cell carcinoma' => 'SCC',
	'biliary tract' => 'biliary',
	'organ' => 'tract|system',
	'urothelial cell' => 'urothelial',
	'leukemia' => 'leukaemia',
	'glioblastoma multiforme' => 'glioblastoma',
	'glioblastoma' => 'glioblastoma multiforme|astrocytoma WHO grade IV|high-grade glioma|high grade glioma',
	'lung small cell carcinoma' => 'small-cell lung carcinoma|small-cell lung cancer|Carcinoma, Small Cell',
	'lung non-small cell carcinoma' => 'non small cell lung cancer|Carcinoma, Non-Small-Cell Lung|non-small-cell lung cancer|non-small cell lung cancer|Non-Small-Cell Lung',
	'hematologic' => 'hematologic|haematologic|hematological|haematological',
	'malignant mesothelioma' => 'mesothelioma',
    'female reproductive organ' => 'gynecologic',
    'adrenocortical carcinoma' => 'Carcinoma, Adrenal Cortical'
);


sub expand_synonyms {
	my $srcfile =  POTTRConfig::get_first_path("data", "disease-ontology-file"); # 'data/doid.obo';
    my $data = join('', file($srcfile));
    
    my @parts = split /\n\n/, $data;
    @parts = grep {/Term/} @parts;
	unshift @parts, @POTTR_parts;
	
	my @expanded_doid_lines;
    for my $p (@parts) {
		my ($id) = ($p =~ /id: (\S+)/);
		my ($name) = ($p =~ /name: (.+)/);
		next if $name =~ /obsolete/;
		next if $p !~ /DO_cancer_slim|TopNodes_DOcancerslim|ICD10CM:C|NCIthesaurus|Solid/; # 
		
		my @synonyms =  ( $p =~ /synonym: "(.+?)".*EXACT/g);
		my @synonyms_new ;
		
	# 	print "\e[1;36m".join("\t", $id, $name, @synonyms, @is_a)."\e[0m\n";
		
		for my $n ($name, @synonyms) {
	# 		$n =~ s///;
			push @synonyms_new, lc($n) if ! exists $DO_id{lc($n)};
			
			# Some dirty hack to include all spelling of names
			for my $i (0..$#synonym_cancer) {
				for my $j (0..$#synonym_cancer) {
					next if $i == $j;
					my $m = lc($n); 
					if ( $m =~ s/(?:^|\b)$synonym_cancer[$i]/$synonym_cancer[$j]/i ) {
						push @synonyms_new, lc($m);
					}
				}
			}
		
			while( my ($key, $values) = each (%synonym_map) ) {
				for my $value ( split /\|/, $values ) {
					my $m = lc($n); 
					if ( $m =~ s/(?:^|\b)$key/$value/i ) {
						push @synonyms_new, lc($m);
					}
				
					for my $i (0..$#synonym_cancer) {
						for my $j (0..$#synonym_cancer) {
							next if $i == $j;
							my $m1 = lc($m); 
							if ( $m1 =~ s/(?:^|\b)$synonym_cancer[$i]/$synonym_cancer[$j]/i ) {
								push @synonyms_new, lc($m1);
							}
						}
					}
				}
			}
		}
		@synonyms_new = uniq(@synonyms_new);
		
		push @expanded_doid_lines, $p."\n";
		push @expanded_doid_lines, map { "synonym: \"$_\" EXACT\n" } @synonyms_new;
		push @expanded_doid_lines, "\n\n";
    }
    
    my $expanded_doid_file = "$srcfile.expanded.obo";
    open FOUT, ">$expanded_doid_file";
    print FOUT @expanded_doid_lines;
    close FOUT;
}


sub cmp_doid {
	my ($x, $y) = @_;
	print "\e[1;33m".join("\t", $x, $y)."\e[0m\n" if $x =~ /POTTR/ and $y =~ /DOID/;
	return 1  if $x =~ /POTTR/ and $y =~ /DOID/;
	return -1 if $x =~ /DOID/ and $y =~ /POTTR/;
# 	my ($xn) = ( $x =~ /:(\d+)/ );
# 	my ($yn) = ( $y =~ /:(\d+)/ );
# 	print "\e[1;33m".join("\t", $x, $y)."\e[0m\n" if $yn < $xn;
	return 0;
# 	return $xn <=> $yn;
}

sub load_DO {
	my $srcfile = POTTRConfig::get_first_path("data", "disease-ontology-file") // 'data/doid.obo'; # default to disease ontology source: data/doid.obo
    my $expanded_doid_file = "$srcfile.expanded.obo";
	
	if ( (! -f $expanded_doid_file) or (mtime( __FILE__ ) > mtime($expanded_doid_file)) ) {
		expand_synonyms;
	}
	
    my $data = join('', file($expanded_doid_file));
    my @parts = split /\n\n+/, $data;
    @parts = grep {/Term/} @parts;
	unshift @parts, @POTTR_parts;
	
    for my $p (@parts) {
	my ($id) = ($p =~ /id: (\S+)/);
	my ($name) = ($p =~ /name: (.+)/);
	next if $name =~ /obsolete/;
	next if $p !~ /DO_cancer_slim|TopNodes_DOcancerslim|ICD10CM:C|NCIthesaurus|Solid/; # 
	
	my @synonyms =  ( $p =~ /synonym: "(.+?)".*EXACT/g);
	my @is_a     =  ( $p =~ /is_a: (\S+)/g);
	
	$DO_is_a{$id}{$_} = 1 for @is_a ;
# 	print "\e[1;36m".join("\t", $id, $name, @synonyms, @is_a)."\e[0m\n";
	$DO_id{lc($name)} = $id;
	for my $n ($name, @synonyms) {
# 		$n =~ s///;
		if ( exists $DO_id{lc($n)} ) {
# 			print "\e[1;32m".join("\t", lc($n))."\e[0m\n" ;
		
			my $old_id = $DO_id{lc($n)};
			$DO_id{lc($n)} = $id if ( cmp_doid($id, $old_id) > 0 ) ;
		} else {
			$DO_id{lc($n)} = $id 
		}
	}
	
	$DO_name{$id} = $name;
# 	print "$p\n";
    }
#    print map { "$_\n----\n" } @parts;    
}


# expand_synonyms;

our $DO_regex ; 

our %DO_cache;

my $DO_cache_file ;


sub get_DO_regex {
	&ON_DEMAND_INIT;
	return $DO_regex;
}

sub DO_match_catype {
	&ON_DEMAND_INIT;
	my $x = shift;
	return @{ $DO_cache{$x} } if ( exists $DO_cache{$x} ) ;
	
	my @matches;
	while ($x =~ /(?:^|\b)\s*($DO_regex)s?\s*(?:\b|$)/ig) {
		my $match = $1;
# 		print join("\t", "\e[1;32m".$DO_id{lc($match)}, "\e[1;33m".$match, "\e[1;36m".$DO_name{ $DO_id{lc($match)} },)."\e[0m\n";	
		push @matches, $DO_id{lc($match)} if exists $DO_id{lc($match)} ;
	}
	@matches = uniq (@matches);
	$DO_cache{$x} = \@matches;
	open  DOCACHE, ">>$DO_cache_file";
	print DOCACHE join("\t", $x, @matches)."\n";
	close DOCACHE;
	chmod 0666, $DO_cache_file;
	return @matches ;
}

sub ngram { return scalar( split /\s+/, $_[0] ) }

sub DO_match_catype_whole_word {
	&ON_DEMAND_INIT;
	my $x = shift;
	
	return $DO_cache{$x}[0] if ( exists $DO_cache{$x} ) ;
	
	my $matched_DOID ;
	while ($x =~ /(?:^|\b)($DO_regex)s?(?:\b|$)/ig) {
		my $match = $1;
		my $n1 = ngram($match);
		my $n2 = ngram($x);
# 		print "\t\t$n1\t$match\n";
# 		print "\t\t$n2\t$x\n";
		
		if ( ($n1 <= $n2 * 0.5) or ( ($n1==$n2) and (length($match) < (length($x) * 0.5)) ) ) { # consider unmatched if just matched < 50% of words # ( ($n1 == $n2) and ( length($match) < (length($x) * 0.5) ) ) or 
			next if $n1 == 1 and $match =~ /^(?:(?:adeno)?carcinoma|cancer|neoplasm|tumou?r|disease)s?$/i; # return $x 
		}
		$matched_DOID = $DO_id{lc($match)} if exists $DO_id{lc($match)} ;
		$DO_cache{$x} = [ ($matched_DOID) ]; #  \@matches;
		open  DOCACHE, ">>$DO_cache_file";
		print DOCACHE join("\t", $x, $matched_DOID)."\n";
		close DOCACHE;
		chmod 0666, $DO_cache_file;
		return $matched_DOID ;
	}
	return $x ;
}

sub get_ancestors {
	&ON_DEMAND_INIT;
	my $x = shift;
	my %is_a = map { $_ => 1 } ( map { /^(?:POTTR|DOID):/ ? ($x) : DO_match_catype($x) } ($x) ); # keys %{ $DO_is_a{$match_id} };
	my @is_a_new ;
	do {@is_a_new = ();
		for my $m (keys %is_a) {
			push @is_a_new, ( grep { ! exists $is_a{$_} } keys %{ $DO_is_a{$m} } ); # exists $DO_name{$_} and 
		}
		$is_a{$_} = 1 for @is_a_new;
	} while (scalar(@is_a_new) > 0);
	my @is_a = sort grep { exists $DO_name{$_} and defined $DO_name{$_} } keys %is_a ;
	return @is_a ;
}

sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;
	
	load_DO ;

	$DO_cache_file = '/tmp/DO_cache.txt';
	chown 0666, $DO_cache_file ;

	if ( -f $DO_cache_file ) {
		for ( file($DO_cache_file) ) {
			chomp;
			next if /^\s*$/;
			my ($x, @matches) = split /\t/, $_;
			$DO_cache{$x} = \@matches;
		}
	}
	
	$DO_regex = join("|", (sort { length($b) <=> length($a) } keys %DO_id) );
}



# print $DO_regex; 

1;
