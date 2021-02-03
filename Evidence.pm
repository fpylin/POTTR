#!/usr/bin/perl
##############################################################################
#
# Evidence.pm - precision oncology trial and therapy recommender
# 
# Evidence integration functions
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

package Evidence;

use strict;
use warnings;

use POSIX;

use lib '.';
use Therapy;
use TSV;

use DOID;
use Rules;

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

sub uniq { my %a; $a{$_} = 1 for(@_) ; return keys %a; }
sub file { open F, "<$_[0]" or die "Unable to open file $_[0].\n"; my @lines = <F>; close F; return @lines; }

our $f_initiailised = undef;

our @debug_msg;
our @error_msg;

# from COSMIC
my @cancer_gene_census = qw( 
A1CF ABI1 ABL1 ABL2 ACKR3 ACSL3 ACSL6 ACVR1 ACVR2A AFDN AFF1 AFF3 AFF4 AKAP9 AKT1 AKT2 AKT3 ALDH2 ALK AMER1 ANK1 APC APOBEC3B 
AR ARAF ARHGAP26 ARHGAP5 ARHGEF10 ARHGEF10L ARHGEF12 ARID1A ARID1B ARID2 ARNT ASPSCR1 ASXL1 ASXL2 ATF1 ATIC ATM ATP1A1 ATP2B3 ATR ATRX AXIN1 AXIN2
B2M BAP1 BARD1 BAX BAZ1A BCL10 BCL11A BCL11B BCL2 BCL2L12 BCL3 BCL6 BCL7A BCL9 BCL9L BCLAF1 BCOR BCORL1 BCR BIRC3 BIRC6
BLM BMP5 BMPR1A BRAF BRCA1 BRCA2 BRD3 BRD4 BRIP1 BTG1 BTK BUB1B
C15orf65 CACNA1D CALR CAMTA1 CANT1 CARD11 CARS CASP3 CASP8 CASP9 CBFA2T3 CBFB CBL CBLB CBLC CCDC6 CCNB1IP1 CCNC CCND1 CCND2 CCND3 CCNE1
CCR4 CCR7 CD209 CD274 CD28 CD74 CD79A CD79B CDC73 CDH1 CDH10 CDH11 CDH17 CDK12 CDK4 CDK6 CDKN1A CDKN1B CDKN2A CDKN2C CDX2 CEBPA CEP89 CHCHD7
CHD2 CHD4 CHEK2 CHIC2 CHST11 CIC CIITA CLIP1 CLP1 CLTC CLTCL1 CNBD1 CNBP CNOT3 CNTNAP2 CNTRL COL1A1 COL2A1 COL3A1 COX6C CPEB3 CREB1 CREB3L1 CREB3L2 CREBBP
CRLF2 CRNKL1 CRTC1 CRTC3 CSF1R CSF3R CSMD3 CTCF CTNNA2 CTNNB1 CTNND1 CTNND2 CUL3 CUX1 CXCR4 CYLD CYP2C8 CYSLTR2
DAXX DCAF12L2 DCC DCTN1 DDB2 DDIT3 DDR2 DDX10 DDX3X DDX5 DDX6 DEK DGCR8 DICER1 DNAJB1 DNM2 DNMT3A DROSHA DUX4L1
EBF1 ECT2L EED EGFR EIF1AX EIF3E EIF4A2 ELF3 ELF4 ELK4 ELL ELN EML4 EP300 EPAS1 EPHA3 EPHA7 EPS15
ERBB2 ERBB3 ERBB4 ERC1 ERCC2 ERCC3 ERCC4 ERCC5 ERG ESR1 ETNK1 ETV1 ETV4 ETV5 ETV6 EWSR1 EXT1 EXT2 EZH2 EZR
FAM131B FAM135B FAM47C FANCA FANCC FANCD2 FANCE FANCF FANCG FANCL FAS FAT1 FAT3 FAT4 FBLN2 FBXO11 FBXW7 FCGR2B FCRL4 FEN1 FES FEV FGFR1 FGFR1OP FGFR2 FGFR3 FGFR4
FH FHIT FIP1L1 FKBP9 FLCN FLI1 FLNA FLT3 FLT4 FNBP1 FOXA1 FOXL2 FOXO1 FOXO3 FOXO4 FOXP1 FOXR1 FSTL3 FUBP1 FUS
GAS7 GATA1 GATA2 GATA3 GLI1 GMPS GNA11 GNAQ GNAS GOLGA5 GOPC GPC3 GPC5 GPHN GRIN2A GRM3
H3F3A H3F3B HERPUD1 HEY1 HIF1A HIP1 HIST1H3B HIST1H4I HLA-A HLF HMGA1 HMGA2 HMGN2P46 HNF1A HNRNPA2B1 HOOK3 HOXA11 
HOXA13 HOXA9 HOXC11 HOXC13 HOXD11 HOXD13 HRAS HSP90AA1 HSP90AB1
ID3 IDH1 IDH2 IGF2BP2 IGH IGK IGL IKBKB IKZF1 IL2 IL21R IL6ST IL7R IRF4 IRS4 ISX ITGAV ITK
JAK1 JAK2 JAK3 JAZF1 JUN
KAT6A KAT6B KAT7 KCNJ5 KDM5A KDM5C KDM6A KDR KDSR KEAP1 KIAA1549 KIF5B KIT KLF4 KLF6 KLK2 KMT2A KMT2C KMT2D KNL1 KNSTRN KRAS KTN1
LARP4B LASP1 LATS1 LATS2 LCK LCP1 LEF1 LEPROTL1 LHFPL6 LIFR LMNA LMO1 LMO2 LPP LRIG3 LRP1B LSM14A LYL1 LZTR1
MACC1 MAF MAFB MALAT1 MALT1 MAML2 MAP2K1 MAP2K2 MAP2K4 MAP3K1 MAP3K13 MAPK1 MAX MB21D2 MDM2 MDM4 MDS2 MECOM MED12 MEN1 MET MGMT MITF MLF1 MLH1
MLLT1 MLLT10 MLLT11 MLLT3 MLLT6 MN1 MNX1 MPL MRTFA MSH2 MSH6 MSI2 MSN MTCP1 MTOR MUC1 MUC16 MUC4 MUTYH MYB MYC MYCL MYCN MYD88 MYH11 MYH9 MYO5A MYOD1
N4BP2 NAB2 NACA NBEA NBN NCKIPSD NCOA1 NCOA2 NCOA4 NCOR1 NCOR2 NDRG1 NF1 NF2 NFATC2 NFE2L2 NFIB NFKB2 NFKBIE NIN
NKX2-1 NONO NOTCH1 NOTCH2 NPM1 NR4A3 NRAS NRG1 NSD1 NSD2 NSD3 NT5C2 NTHL1 NTRK1 NTRK3 NUMA1 NUP214 NUP98 NUTM1 NUTM2B NUTM2D
OLIG2 OMD P2RY8
PABPC1 PAFAH1B2 PALB2 PATZ1 PAX3 PAX5 PAX7 PAX8 PBRM1 PBX1 PCBP1 PCM1 PDCD1LG2 PDE4DIP PDGFB PDGFRA PDGFRB PER1 PHF6 PHOX2B PICALM PIK3CA PIK3CB PIK3R1 PIM1 PLAG1
PLCG1 PML PMS1 PMS2 POLD1 POLE POLG POLQ POT1 POU2AF1 POU5F1 PPARG PPFIBP1 PPM1D PPP2R1A PPP6C PRCC PRDM1 PRDM16 PRDM2 PREX2 PRF1 PRKACA PRKAR1A PRKCB
PRPF40B PRRX1 PSIP1 PTCH1 PTEN PTK6 PTPN11 PTPN13 PTPN6 PTPRB PTPRC PTPRD PTPRK PTPRT PWWP2A
QKI
RABEP1 RAC1 RAD17 RAD21 RAD50 RAD51B RAF1 RALGDS RANBP2 RAP1GDS1 RARA RB1 RBM10 RBM15 RECQL4 REL RET RFWD3 RGPD3 RGS7 RHOA RHOH RMI2 RNF213 RNF43 
ROBO2 ROS1 RPL10 RPL22 RPL5 RPN1 RSPO2 RSPO3 RUNX1 RUNX1T1
S100A7 SALL4 SBDS SDC4 SDHA SDHAF2 SDHB SDHC SDHD SEPT5 SEPT6 SEPT9 SET SETBP1 SETD1B SETD2 SETDB1 SF3B1 SFPQ SFRP4 SGK1 SH2B3 SH3GL1 SHTN1 SIRPA SIX1 SIX2 SKI SLC34A2 SLC45A3
SMAD2 SMAD3 SMAD4 SMARCA4 SMARCB1 SMARCD1 SMARCE1 SMC1A SMO SND1 SNX29 SOCS1 SOX2 SOX21 SPECC1 SPEN SPOP SRC SRGAP3 SRSF2 SRSF3
SS18 SS18L1 SSX1 SSX2 SSX4 STAG1 STAG2 STAT3 STAT5B STAT6 STIL STK11 STRN SUFU SUZ12 SYK
TAF15 TAL1 TAL2 TBL1XR1 TBX3 TCEA1 TCF12 TCF3 TCF7L2 TCL1A TEC TENT5C TERT TET1 TET2 TFE3 TFEB TFG TFPT TFRC TGFBR2 THRAP3 TLX1 TLX3 TMEM127 TMPRSS2 
TNC TNFAIP3 TNFRSF14 TNFRSF17 TOP1 TP53 TP63 TPM3 TPM4 TPR TRA TRAF7 TRB TRD TRIM24 TRIM27 TRIM33 TRIP11 TRRAP TSC1 TSC2 TSHR
U2AF1 UBR5 USP44 USP6 USP8
VAV1 VHL VTI1A
WAS WDCP WIF1 WNK2 WRN WT1 WWTR1
XPA XPC XPO1
YWHAE
ZBTB16 ZCCHC8 ZEB1 ZFHX3 ZMYM2 ZMYM3 ZNF331 ZNF384 ZNF429 ZNF479 ZNF521 ZNRF3 ZRSR2
);

my @oncogene_list = qw(
A1CF ABL1 ABL2 ACKR3 ACVR1 AFDN AFF3 AFF4 AKT1 AKT2 AKT3 ALK APOBEC3B AR ARAF ARHGAP5 ARNT ATF1 ATP1A1 
BCL11A BCL11B BCL2 BCL2L12 BCL3 BCL6 BCL9 BCL9L BCORL1 BIRC3 BIRC6 BMPR1A BRAF BRD3 BRD4 BTK 
CACNA1D CALR CARD11 CBL CBLC CCND1 CCND2 CCND3 CCNE1 CCR4 CCR7 CD28 CD74 CD79A CD79B CDH17 CDK4 CDK6 CDKN1A CHD4 CHST11 CIC CREB1 CREB3L2 CREBBP CRLF2 CRTC1 CSF1R CSF3R CTNNA2 CTNNB1 CTNND2 CUX1 CXCR4 CYSLTR2 
DAXX DDB2 DDIT3 DDR2 DDX5 DDX6 DEK DGCR8 
EGFR ELF4 ELK4 EPAS1 ERBB2 ERBB3 ERBB4 ERG ESR1 ETV1 ETV4 ETV5 EWSR1 EZH2 
FCGR2B FCRL4 FES FEV FGFR1 FGFR2 FGFR3 FGFR4 FLI1 FLT3 FLT4 FOXA1 FOXL2 FOXO1 FOXO3 FOXO4 FOXP1 FOXR1 FSTL3 FUBP1 
GATA1 GATA2 GATA3 GLI1 GNA11 GNAQ GNAS GPC3 GRM3 
H3F3A H3F3B HEY1 HIF1A HIP1 HIST1H3B HLF HMGA1 HMGA2 HNRNPA2B1 HOXA11 HOXA13 HOXA9 HOXC11 HOXC13 HOXD11 HOXD13 HRAS 
IDH1 IDH2 IKBKB IL6ST IL7R IRF4 IRS4 
JAK1 JAK2 JAK3 JUN 
KAT6A KAT7 KCNJ5 KDM5A KDM6A KDR KIT KLF4 KMT2A KMT2D KNSTRN KRAS 
LCK LEF1 LMO1 LMO2 LPP LYL1 
MACC1 MAF MAFB MALAT1 MALT1 MAML2 MAP2K1 MAP2K2 MAP2K4 MAP3K1 MAP3K13 MAPK1 MDM2 MDM4 MECOM MET MITF MLLT10 MN1 MPL MRTFA MSI2 MTCP1 MTOR MUC16 MUC4 MYB MYC MYCL MYCN MYD88 MYOD1 NCOA2 
NFATC2 NFE2L2 NFKB2 NKX2-1 NOTCH1 NOTCH2 NPM1 NR4A3 NRAS NSD2 NSD3 NT5C2 NTRK1 NTRK3 NUP98 NUTM1 
OLIG2 
P2RY8 PABPC1 PAX3 PAX5 PBX1 PDCD1LG2 PDGFB PDGFRA PDGFRB PIK3CA PIK3CB PIM1 PLAG1 PLCG1 POLQ POU2AF1 POU5F1 PPM1D PRDM16 PREX2 PRKACA PRKAR1A PSIP1 PTK6 PTPN11 
QKI 
RAC1 RAD21 RAF1 RAP1GDS1 RARA RECQL4 REL RET RHOA ROS1 RSPO3 RUNX1 RUNX1T1 
SALL4 SET SETBP1 SETDB1 SF3B1 SGK1 SH3GL1 SIX1 SIX2 SKI SMO SND1 SOX2 SRC SRSF2 SRSF3 SSX1 SSX2 SSX4 STAT3 STAT5B STAT6 STIL SUZ12 SYK 
TAF15 TAL1 TAL2 TBL1XR1 TBX3 TCF3 TCF7L2 TCL1A TEC TERT TET1 TFE3 TFEB TLX1 TLX3 TNC TNFRSF17 TP63 TRIM24 TRIM27 TRRAP TSHR 
U2AF1 UBR5 USP6 USP8 
WAS WT1 WWTR1 
XPO1 
ZEB1 ZNF521
);

my @TSG_list = qw(
ABI1 ACVR2A AMER1 APC APOBEC3B ARHGAP26 ARHGEF10 ARHGEF10L ARHGEF12 ARID1A ARID1B ARID2 ARNT ASXL1 ASXL2 ATM ATP1A1 ATP2B3 ATR ATRX AXIN1 AXIN2 
B2M BAP1 BARD1 BAX BAZ1A BCL10 BCL11B BCL9L BCOR BCORL1 BIRC3 BLM BMPR1A BRCA1 BRCA2 BRIP1 BTG1 BTK BUB1B 
CAMTA1 CARS CASP3 CASP8 CASP9 CBFA2T3 CBFB CBL CBLB CBLC CCDC6 CCNB1IP1 CCNC CD274 CDC73 CDH1 CDH10 CDH11 CDK12 CDKN1A CDKN1B CDKN2A CDKN2C CDX2 CEBPA CHD2 CHEK2 CIC CIITA CLTC CLTCL1 CNBP CNOT3 CNTNAP2 CPEB3 CREB3L1 CREBBP CSMD3 CTCF CUL3 CUX1 CYLD 
DAXX DDB2 DDX10 DDX3X DICER1 DNM2 DNMT3A DROSHA 
EBF1 EED EIF3E ELF3 ELF4 ELL EP300 EPAS1 EPS15 ERBB4 ERCC2 ERCC3 ERCC4 ERCC5 ESR1 ETNK1 ETV6 EXT1 EXT2 EZH2 
FANCA FANCC FANCD2 FANCE FANCF FANCG FANCL FAS FAT1 FAT4 FBLN2 FBXO11 FBXW7 FEN1 FES FH FHIT FLCN FOXL2 FOXO1 FOXO3 FOXO4 FUS 
GATA1 GATA3 GPC3 GPC5 GRIN2A 
HNF1A HOXA11 HOXA9 
ID3 IGF2BP2 IKZF1 IRF4 IRS4 JAK1 
KAT6B KDM5C KDM6A KEAP1 KLF4 KLF6 KMT2C KMT2D KNL1 
LARP4B LATS1 LATS2 LEF1 LEPROTL1 LRIG3 LRP1B LZTR1 
MALAT1 MAP2K4 MAP3K1 MAP3K13 MAX MED12 MEN1 MGMT MLF1 MLH1 MRTFA MSH2 MSH6 MUTYH MYH9 
N4BP2 NAB2 NBN NCOA4 NCOR1 NCOR2 NDRG1 NF1 NF2 NFE2L2 NFKB2 NFKBIE NKX2-1 NOTCH1 NOTCH2 NRG1 NTHL1 NTRK1 
PABPC1 PALB2 PATZ1 PAX5 PBRM1 PER1 PHF6 PHOX2B PIK3R1 PML PMS2 POLD1 POLE POLG POLQ POT1 PPARG PPP2R1A PPP6C PRDM1 PRDM2 PRF1 PRKAR1A PTCH1 PTEN PTK6 PTPN13 PTPN6 PTPRB PTPRC PTPRD PTPRK PTPRT 
QKI 
RAD17 RAD21 RAD50 RAD51B RANBP2 RB1 RBM10 RECQL4 RFWD3 RHOA RHOH RMI2 RNF43 ROBO2 RPL10 RPL22 RPL5 RSPO2 RUNX1 RUNX1T1 
SBDS SDHA SDHAF2 SDHB SDHC SDHD SETD1B SETD2 SFPQ SFRP4 SH2B3 SIRPA SLC34A2 SMAD2 SMAD3 SMAD4 SMARCA4 SMARCB1 SMARCD1 SMARCE1 SMC1A SOCS1 SOX21 SPEN SPOP STAG1 STAG2 STAT5B STK11 SUFU SUZ12 
TBL1XR1 TBX3 TCF3 TENT5C TERT TET1 TET2 TGFBR2 TMEM127 TNFAIP3 TNFRSF14 TP53 TP63 TPM3 TRAF7 TRIM24 TRIM33 TSC1 TSC2 
USP44 
VHL 
WIF1 WNK2 WRN WT1 
XPA XPC 
YWHAE 
ZBTB16 ZFHX3 ZMYM3 ZNF331 ZNRF3 ZRSR2
);

my %cancer_gene_census = map { $_ => 1 } @cancer_gene_census;

my %other_biomarkers = (
	'AURKA'   => 'TSG',
	'BCL2L1'  => 'oncogene',    'BCL2L11'=> 'TSG',      'BRD2'   => 'TSG',       'BRDT'   => 'TSG', 
	'CD20'    => 'Other',       'CHEK1'  => 'TSG',      'CEACAM5' => 'Other',    'CD33'   => 'Other',  'CD38' => 'Other',  'CD123' => 'Other',  'CD19' => 'Other', 'CD30' => 'Other',  'CXCL13' => 'Other',
	'EMSY'    => 'oncogene',    'ERRFI1'  => 'TSG',     'ERCC1'   => 'TSG', 
	'FANCI'   => 'TSG',         'FANCM'   => 'TSG',     'FGF19'   => 'oncogene', 'FOLR' => 'Other', 'FRS2' => 'oncogene',
	'GLI2'    => 'TSG',
	'HGF'     => 'oncogene',    
	'MCL1'    => 'oncogene',   'MSLN'    => 'protein',  'MTAP'   => 'TSG',  
	'NECTIN4' => 'Other',      'NOTCH3'  => 'oncogene', 'NTRK2'  => 'oncogene', 'NTRK3'  => 'oncogene',
	'PDGFB'   => 'oncogene',   'PPP2R2A' => 'TSG',      'PRKCA'  => 'oncogene', 'PRKC'   => 'oncogene', 'PRMT1'  => 'TSG',      'PSMA'    => 'protein', 'PGR' => 'oncogene', 
	'RAD51C'  => 'TSG',        'RAD51D'  => 'TSG',      'RAD54L' => 'TSG',      'RICTOR' => 'oncogene',
	'SMARCA1' => 'TSG',        'SLAMF7' => 'Other',     'SSTR2'  => 'Other',    'SMARCA2' => 'TSG',     'SLC1A5' => 'Other',
	'TACSTD2' => 'Other',
);

our %cancer_gene_type ;

$cancer_gene_type{$_} = 'TSG' for @TSG_list;
$cancer_gene_type{$_} = 'oncogene' for @oncogene_list;
$cancer_gene_type{ $other_biomarkers{$_} } = $other_biomarkers{$_} for keys %other_biomarkers;



sub deutf { my $x = shift; $x =~ s/\xc2\xa0//g; return $x ; }
sub hl { my ($col, $x) = @_; return "\e[1;${col}m$x\e[0m"; }

sub encode_alteration_proper {
	my $x = shift;
	my $bm = '';
	$bm = $1 if ( $x =~ s/^([A-Z\-a-z0-9 ]+:)// );
# 	print "\e[1;42;35m".$bm.'|'.$x."\e[0m\n";

	if ( $bm =~ /^((?:Microsatellite.Instability|Tumour.Mutation(?:al)?.Burden|(?:Loss-of-heterozygosity|Homologous Recombination Deficiency)|Homologous.Recombination.Deficiency|Mismatch.repair):)$/i ) { 
		$bm = lc($1);
		$bm =~ s/ /_/g;
	}


	$x =~ s/mut$//;
	$x =~ s/^\s+|\s+$//g;
	if ( $x =~ /^(amplification|overexpression|loss of (?:protein )?expression|(?:homozygous )deletion|wildtype|(?:oncogenic|truncating) mutation|fusion|internal tandem duplication|kinase domain duplication|(?:loss|gain)-of-function mutation|.*variant|alteration|high|deficient)s?$/i ) { 
		$x = lc($1) ;
	} elsif ( $x =~ /^((?:DNA binding|kinase) domain (?:deletion|insertion|duplication|(?:missense )?mutation))s?$/i ) { 
		$x = lc($1);
	} elsif ( $x =~ /^(Exon \d+ (?i:deletion|(?:splic\w+ |skipping )?mutation|insertion|insertions\/deletions|indel))s?$/i ) { 
		$x = lc($1);
		$x =~ s/insertions\/deletions/indel/i;
	} elsif ( $x =~ /^((?:protein )?(?:over)?expression)s?$/i ) { 
		$x = lc($1);
	} elsif ( $x =~ /^(.*) (fusion)$/i ) { 
		$x = $1."_".lc($2);
	} elsif ( $x =~ /^[A-Z](\d+)ins$/i ) { 
		$x = lc("codon_$1_inframe_insertion");
	} else {
		return ($bm // '').$x;
	}
	$x =~ s/ /_/g;
	return ($bm // '').$x;
}

sub encode_alteration {
	my ($biomarker, $alteration) = @_;
	
	my @parts = split /(?: and )/i, $alteration;
	my @retval;
	
	for my $x (@parts) {
		while ( $biomarker =~ /(Microsatellite.Instability|Tumour.Mutation(?:al)?.Burden|(?:Loss-of-heterozygosity|Homologous Recombination Deficiency)|Homologous.Recombination.Deficiency|Mismatch.repair)/ig ) {
			my $m = $1;
			if ( $x =~ /$m:/i ) {
				$x = lc($x);
			}
		}

		my $f_neg = 0;
		my $f_germline = 0;
		$f_neg = 1 if  $x =~ s/^NOT +//gi ;  
		$f_germline = 1 if  $x =~ s/\s*\(\s*germline\s*\)//g ;  
		
		if ( $x =~ /\s+\+\s+/ ) {
			my @xparts = split /\s+\+\s+/, $x ;
			push @retval, join('; ', ( map { ($f_neg?'NOT ':'').encode_alteration_proper($_) } @xparts ) );
			next;
		}
		
		push @retval, ($f_neg?'NOT ':'').encode_alteration_proper($x).($f_germline ? ',germline' : '');
	}
	
	my $f_biomarker_has_plus = ( $biomarker =~ /\+/ );
	
	return join("; ", ( map { 
			( /^((?i:not)) +(?<bm>.*)/ ? 
			  ( "$1 ".($f_biomarker_has_plus  ? "" : "$biomarker:").$+{bm} ) 
			  :
			  ( ( $f_biomarker_has_plus ? "" : "$biomarker:"). $_ ) 
			)
			} @retval 
		) 
	);
}

sub expand_alterations {
	my $x = shift;
	my @x = ($x);
# 	push @x, 'fusion' if $x =~ / fusion/i;
	return @x;
}

sub process_tumour_type {
	my $x = shift;
	$x =~ s/Glioma/Glioblastoma/;
	return $x;
}


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

sub check_biomarker {
	my $biomarker = shift;
	for ( split /\s*\+\s*/, $biomarker ) {
		return undef if /Microsatellite Instability|Tumour Mutation(?:al)? Burden|(?:Homologous Recombination Deficiency|Loss-of-heterozygosity) Score|Mismatch repair/i;
		push @error_msg, hl(31, "ERROR: Biomarker ''$biomarker'' not found\n") if ! exists $cancer_gene_census{$_} and ! exists $other_biomarkers{$_};
	}
}

our %data_d;
our %data_c;

sub mtime { my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($_[0]); return $mtime; }

sub gen_rule_knowledge_base {
	&ON_DEMAND_INIT;
	my $srckb = shift; # input file name
	my $srcf = shift;  # source file
	my $outf = shift;  # cache file
	
	my @all_rules_treatments;
	
	my $TSV_master = TSV->new($srcf);
	
	%data_d = ();
	%data_c = ();

	my @fields = @{ $TSV_master->{'fields'} };

	my ($fname_biomarker)   = grep { /Biomarker|Hugo\s+Symbol/x } @fields ;
	my ($fname_alteration)  = grep { /Alteration|Alteration/x } @fields ;
	my ($fname_cancer_type) = grep { /Tumour\s+Type|Cancer\s+Type/x } @fields ;
	my ($fname_evidence)    = grep { /Evidence|PMIDs\s+for\s+drug/x } @fields ;
	my ($fname_drugs)       = grep { /Drugs|Drug\(s\)/x } @fields ;
	my ($fname_tier)        = grep { /Tier|Level/x } @fields ;
	my $mtime_db            = strftime("%d/%m/%Y", localtime( mtime($srcf) ) );

	# First, register all known drug combinations from the database
	for my $row ( @{ $TSV_master->{'data'} } ) {
		my $treatment = deutf( $$row{$fname_drugs} );
		my @treatments = split /\s*(?:[;]|,)\s+/, $treatment;
		for my $rx ( @treatments ) {
			if ( Therapy::is_combination( $rx ) ) {
				Therapy::register_known_combination( $rx );
			}
		}
	}

	# Load repurposing rules:
	my %repurposing_retier_cancer_type ;
	my %repurposing_retier_related_mutation ;
	
	my ($repurposing_retier_cancer_type_str) = POTTRConfig::get('repurposing-retier:cancer-type');
	if ( defined $repurposing_retier_cancer_type_str ) {
		push @debug_msg, "repurposing-retier:cancer-type = $repurposing_retier_cancer_type_str\n";
		%repurposing_retier_cancer_type = map { my ($a, $b) = split /\s*:\s*/, $_; $a => $b } split /\s*,\s*/, $repurposing_retier_cancer_type_str;
	}

	my ($repurposing_retier_related_mutation_str) = POTTRConfig::get('repurposing-retier:related-mutation');
	if ( defined $repurposing_retier_related_mutation_str ) {
		push @debug_msg, "repurposing-retier:related-mutation = $repurposing_retier_related_mutation_str\n";
		%repurposing_retier_related_mutation = map { my ($a, $b) = split /\s*:\s*/, $_; $a => $b } split /\s*,\s*/, $repurposing_retier_related_mutation_str;
	}

	# Then, check if there are rules already cached and uptodate. If so, use it.
	if ( defined($outf) and ( -f $outf ) and ( mtime($srcf) < mtime($outf) ) and ( mtime(__FILE__) < mtime($outf) ) ) { 
		return file($outf);
	}
	
	# Otherwise generate new rules.
	
	
	my %rules_treatment_by_catype_mutation_position;
	
	my $line = 1;
	for my $row ( @{ $TSV_master->{'data'} } ) {
		$line ++;
		my $source      = $$row{'Source'} // '';   $source  =~ s/-.*//;
		my $date        = $$row{'Date'} // $mtime_db;
		my $tier        = $$row{$fname_tier};
		my $biomarker   = deutf( $$row{$fname_biomarker} );             $biomarker =~ s/^\s*|\s*$//;
		my $alterations = deutf( $$row{$fname_alteration} );
		my $tumour_type = deutf( $$row{$fname_cancer_type} );
		my $treatment   = deutf( $$row{$fname_drugs} );
		my $evidence    = deutf( $$row{$fname_evidence} // '' );        $evidence =~ s/;/,/g;
		my $comments    = deutf( $$row{'Comments'} // '' );
		my $treatment_class = Therapy::get_treatment_class( $$row{$fname_drugs}, $biomarker );
		
		my $version_str = ($date =~ m|(\d+)/(\d+)/(\d+)| ? "KBver:$3$2$1" : '');
		
		my $tier_partial_match_catype = $repurposing_retier_cancer_type{$tier} // 'U';
		
		push @debug_msg, "\e[1;41;37mUntranslated tier $tier\e[0m\n" if $tier_partial_match_catype eq 'U' ;
		
# 		( ( $tier =~ /^R/) ? 'R2B' : ( ($tier =~ /^[123]/) ? '3B' : '4B' ) );

		$tumour_type = process_tumour_type($tumour_type);
		
		check_biomarker($biomarker);

		push @debug_msg, "\e[1;37m".join("\t", $source, $date, $tier, $biomarker, $alterations, $tumour_type, $comments, $treatment, $evidence )."\e[0m\n";
		
		my @rules;

		if ( $biomarker =~ /\s/ and $biomarker !~ /\+/ ) {
			$biomarker = lc($biomarker);
			$biomarker =~ s/\s+/_/g ;
		}
		
		
		my @lhs_catypes_neg ;
		if ( $tumour_type =~ /(?:\s*,?\s+)?except\s+/i ) {
			$tumour_type = $`;
# 			print "[[$tumour_type]]\n";
			@lhs_catypes_neg = split /\s*[,;]\s*/, $';
		}
		
		my @lhs_catypes = split /\s*[;]\s*/, $tumour_type;
		my @lhs_catypes_ancestors; # determining cancer type "ancestors" to avoid inference across boundaries
		
		for my $c (@lhs_catypes) {
			my ($catype_match) = DO_match_catype($c);
# 			print ">>\t$tumour_type\t$c\t$catype_match\n";
			my $x = "catype:".(defined $catype_match ? $catype_match : $c) ;
			$c = $x ;
			
			if ( defined $catype_match ) {
				my @ancerstors = DOID::get_ancestors( $catype_match );
				my @lhs_ancestors ;
				for my $a (@ancerstors) {
					push @lhs_ancestors , ( "catype:".$a, "catype:".$DO_name{$a}, "catype_name:$DO_name{$a}" ) ;
# 					print "[$a]\t$DO_name{$a}\n";
				}
				push @lhs_catypes_ancestors, @lhs_ancestors ;
			}
		}
		
		if ( scalar @lhs_catypes_ancestors ) { # FIXME: temporary fix for unrooted cancer type classes, until a cancer type ontology is implemented
			my ($a) = DO_match_catype('Solid tumour');
			($a) = DO_match_catype('Liquid cancer') if ( grep { /ha?ematologic|myelodysplas|macroglobulina?emia|leuka?emia|myeloma/i } @lhs_catypes_ancestors );
# 			print map {"|| $_\n"} @lhs_catypes_ancestors ;
# 			print "\n";
			die $a if ! defined $DO_name{$a};
			push @lhs_catypes_ancestors, ( "catype:".$a, "catype:".$DO_name{$a}, "catype_name:$DO_name{$a}" ) ;
		}
		
		my $lhs_catype = ( (scalar(@lhs_catypes) > 1) ? "(".join(" OR ", @lhs_catypes).")" : $lhs_catypes[0] );
		my $lhs_catype_ancestors = ( (scalar(@lhs_catypes_ancestors) > 1) ? "(".join(" OR ", @lhs_catypes_ancestors).")" : $lhs_catypes_ancestors[0] );
		my $lhs_catype_neg ;
		
		if ( scalar @lhs_catypes_neg ) {
			for my $c (@lhs_catypes_neg) {
				my ($catype_match) = DO_match_catype($c);
				my $cstr = defined $catype_match ? join('; ', ( map { "NOT $_" } ( "catype:".$catype_match, "catype:".$DO_name{$catype_match}, "catype_name:$DO_name{$catype_match}" ) ) )  : "NOT catype:$catype_match";
				$lhs_catype_neg = defined $lhs_catype_neg ? join("; ", $lhs_catype_neg, $cstr) : $cstr ;
			}
		}

		
		my $alterations_neg = '';
		my $alterations_sep = ',;';
		
# 		$alterations_sep = ';';  # if exception is specified, the entire list after exception is interpreted as negative
		if ( $alterations =~ /\s*(?:except|(?:but) not)\s*/i ) {
			($alterations, $alterations_neg) = ($`, $');
		}
		
		
		
		my @alterations = grep { length($_) } map { s/^\s*|\s*$//; $_ } split /\s*[$alterations_sep]\s*/, $alterations ;
		
		my @alterations_neg = grep { length($_) } map { s/^\s*|\s*$//; $_ } ( split /\s*[$alterations_sep]\s*/, $alterations_neg );
		
		@alterations = map {
				join( "; ", encode_alteration($biomarker, $_), ( map { "NOT ".encode_alteration($biomarker, $_) } map { expand_alterations($_) } @alterations_neg ) )
			} ( map { expand_alterations($_) } @alterations )
		;
		
		push @debug_msg, "\e[1;41;37m".join(" ", @alterations )."\e[0m\n" if scalar @alterations_neg ;
		
		for my $alt (@alterations) {
			my $lhs_alteration = $alt;
			my $lhs_alteration_pos ;
			
			$lhs_alteration =~ /^(.+?:)[A-Z]([0-9]+)$/ and do { $lhs_alteration = $lhs_alteration_pos = $1."codon_".$2."_missense_variant"; };
			
			$lhs_alteration =~ /^(.+?:)[A-Z]([0-9]+)[A-Z]$/ and do { 
				$lhs_alteration_pos = $1."codon_".$2."_missense_variant"; 
				$lhs_alteration_pos = undef if $lhs_alteration_pos eq $lhs_alteration ;
			};
			
			my $f_lhs_alteration_is_neg = ( scalar(@alterations) == scalar(grep { /^NOT / } @alterations) );
			
			$lhs_alteration .= '; NOT no_negative_predictive_biomarkers' if $f_lhs_alteration_is_neg  ; 
		
	# 		my $rhs_catype = $tumour_type; # "catype:".$tumour_type;
			
			my @treatments = split /\s*(?:[;]|,)\s+/, $treatment;
			my @rhs_catypes = split /\s*(?:[;]|,)\s+/, $tumour_type;
			
			for my $rhs_catype (@rhs_catypes) {
				for my $rx (@treatments) {
					my $rx = Therapy::get_normalised_treatment_name( $rx );
					my $drug_class_regimen = Therapy::get_treatment_class($rx, $biomarker);
					
					push @debug_msg, join("\t", 
							hl(31, $tier), hl(33, $alt), hl(35, $rx), hl(32, $drug_class_regimen), hl(34,$rhs_catype), hl(36, (my $normalised_class = Therapy::get_normalised_treatment_class_name($rx)) )
						)."\n\n"; 
					
					my $f_no_specific_drug = ( $drug_class_regimen =~ /UNDEF(?:\s+\+\s+UNDEF)*/ );
					
					push @error_msg, "\e[1;31mERROR: Line $line: ''$rx'' not exist in the database\e[0m" if $f_no_specific_drug and ! length $normalised_class ;
					
					$drug_class_regimen = $normalised_class if $f_no_specific_drug ;
					
					my @tags = ( $version_str );
					push @tags, "evidence:$evidence" if length $evidence;
					
					if ( $drug_class_regimen ) {
						my $rhs_treatment_class_str = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $tier, histotype agnostic)";
						
						my $ppalt = $alt; 
						$ppalt =~ s/;/ and /g;
						
						if ( $tumour_type =~ /(?:All)?(?:.*Solid).*Tumou?rs/i ) {
							my $rhs_treatment_class_str = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $tier, histotype agnostic)";
							
							push @rules, mkrule( [$lhs_alteration, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str, "CERTAIN:treatment_class", @tags ) ] );
							
							if ( defined $lhs_alteration_pos ) {
								my $rhs_treatment_class_str = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $repurposing_retier_related_mutation{$tier}, from $ppalt)";
								
								push @{ $rules_treatment_by_catype_mutation_position{$drug_class_regimen}{$lhs_alteration_pos}{$lhs_alteration} }, 
									mkrule( [$lhs_alteration_pos, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str, "INFERRED:treatment_class", @tags ) ] );
							}
							
						} else {
							my $rhs_treatment_class_str = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $tier)";
							
							my $rhs_treatment_class_str_type_specific = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $tier_partial_match_catype, inferred from $rhs_catype)";

							push @rules, mkrule( [$lhs_alteration, $lhs_catype, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str, "CERTAIN:treatment_class", @tags ) ] );
							
							if ( ! $f_lhs_alteration_is_neg ) {
								push @rules, mkrule( [$lhs_alteration, "NOT ".$lhs_catype, $lhs_catype_ancestors, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str_type_specific, "INFERRED:treatment_class", @tags ) ] );
							}
							
							if ( defined $lhs_alteration_pos ) {
								my $rhs_treatment_class_str = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $repurposing_retier_related_mutation{$tier}, from $ppalt)";
								
								my $rhs_treatment_class_str_type_specific = "$biomarker:treatment_class:$drug_class_regimen ($srckb LOE: $repurposing_retier_related_mutation{$tier_partial_match_catype}, inferred from $rhs_catype and $ppalt)";
								
								push @error_msg, "\e[1;31mUndefined tier $tier_partial_match_catype\e[0m\n" if ! exists $repurposing_retier_related_mutation{$tier_partial_match_catype};
							
								push @{ $rules_treatment_by_catype_mutation_position{$drug_class_regimen}{$lhs_alteration_pos}{$lhs_alteration} }, 
									mkrule( [$lhs_alteration_pos, $lhs_catype, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str, "INFERRED:treatment_class", @tags ) ] );
							
								if ( ! $f_lhs_alteration_is_neg ) {
									push @{ $rules_treatment_by_catype_mutation_position{$drug_class_regimen}{$lhs_alteration_pos}{$lhs_alteration} }, 
										mkrule( [$lhs_alteration_pos, "NOT ".$lhs_catype, $lhs_catype_ancestors, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_class_str_type_specific, "INFERRED:treatment_class", @tags ) ] );
								}
							}
						}
						
						# AUDIT:
						for my $alt (@alterations) {
							my $z = $data_c{lc($rhs_catype)}{$biomarker}{$alt}{$drug_class_regimen};
# 							$z = $tier if ! defined $z;
# 							$z = $tier if POTTR::score_tier($tier) < POTTR::score_tier($z) ;
							$data_c{lc($rhs_catype)}{$biomarker}{$alt}{$drug_class_regimen}{$tier}++; #  = $z;
						}
					}


					if ( ! $f_no_specific_drug ) {  # Making rules 
						my $rhs_treatment_str = "$biomarker:treatment:$rx ($srckb LOE: $tier, histotype agnostic)";
						
						my $ppalt = $alt; 
						$ppalt =~ s/;/ and /g;
							
						if ( $tumour_type =~ /(?:All)?(?:.*Solid).*Tumou?rs/i ) {
							my $rhs_treatment_str = "$biomarker:treatment:$rx ($srckb LOE: $tier, histotype agnostic)";
						
							
							push @rules, mkrule( [$lhs_alteration, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str, "CERTAIN:treatment", @tags ) ] );
							
							if ( defined $lhs_alteration_pos ) {
								my $rhs_treatment_str = "$biomarker:treatment:$rx ($srckb LOE: $repurposing_retier_related_mutation{$tier}, histotype agnostic and from $ppalt )";
								
								push @{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos}{$lhs_alteration} }, 
									mkrule( [$lhs_alteration_pos, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str, "INFERRED:treatment", @tags ) ] );
							}
							
						} else {
							my $rhs_treatment_str = "$biomarker:treatment:$rx ($srckb LOE: $tier)";
							
							my $rhs_treatment_str_type_specific = "$biomarker:treatment:$rx ($srckb LOE: $tier_partial_match_catype, inferred from $rhs_catype)";

							push @rules, mkrule( [$lhs_alteration, $lhs_catype, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str, "CERTAIN:treatment", @tags ) ] );
							
							if ( ! $f_lhs_alteration_is_neg ) {
								push @rules, mkrule( [$lhs_alteration, "NOT ".$lhs_catype, $lhs_catype_ancestors, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str_type_specific, "INFERRED:treatment", @tags ) ] );
							}
							
							if ( defined $lhs_alteration_pos ) {
								my $rhs_treatment_str = "$biomarker:treatment:$rx ($srckb LOE: $repurposing_retier_related_mutation{$tier}, from $ppalt )";
								
								my $rhs_treatment_str_type_specific = "$biomarker:treatment:$rx ($srckb LOE: $repurposing_retier_related_mutation{$tier_partial_match_catype}, inferred from $rhs_catype and $ppalt)";
								
								push @{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos}{$lhs_alteration} }, 
									mkrule( [$lhs_alteration_pos, $lhs_catype, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str, "INFERRED:treatment", @tags ) ] );
									
								if ( ! $f_lhs_alteration_is_neg ) {
									push @{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos}{$lhs_alteration} }, 
										mkrule( [$lhs_alteration_pos, "NOT ".$lhs_catype, $lhs_catype_ancestors, $lhs_catype_neg], [ Facts::mk_fact_str( $rhs_treatment_str_type_specific, "INFERRED:treatment", @tags ) ] );
								}
							}
						}
						# AUDIT:
						for my $alt (@alterations) {
							$data_d{lc($rhs_catype)}{$biomarker}{$alt}{$rx}{$tier}++;
						}
					}

					push @debug_msg, map { "$_\n" } @rules;
					push @debug_msg, "\n";
					push @all_rules_treatments, @rules;
				}
			}
		}
		push @debug_msg, "\n";
	}
	
	for my $rx ( sort keys %rules_treatment_by_catype_mutation_position) {
		for my $lhs_alteration_pos ( sort keys %{ $rules_treatment_by_catype_mutation_position{$rx} } ) {
			my $n_alterations = scalar keys %{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos} };
			next if $n_alterations <= 1;
			
			for my $lhs_alteration ( sort keys %{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos} } ) {
				push @debug_msg, "$n_alterations alterations: \t$lhs_alteration_pos\t$lhs_alteration\n";
				my @rules = @{ $rules_treatment_by_catype_mutation_position{$rx}{$lhs_alteration_pos}{$lhs_alteration} };
				@rules = uniq @rules ;
				@rules = sort @rules ;
				push @all_rules_treatments, @rules ;
				for my $rule ( @rules ) {
					push @debug_msg, "$rule\n";
				}
			}
		}
	}
	
	write_array_uniq($outf, @all_rules_treatments) if defined $outf;
	
	return @all_rules_treatments;
}

sub ON_DEMAND_INIT {
	return if $f_initiailised ;
	$f_initiailised = 1;
}

# print @debug_msg;
# print @error_msg;


1;
