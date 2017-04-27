#!/usr/bin/env perl
#===============================================================================
#
#         FILE: Chorus_perl.pl
#        USAGE: ./Chorus_perl.pl
#  DESCRIPTION: design probe using Chorus method
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 04/27/2017 03:24:40 PM
#===============================================================================

use strict;
use warnings;

use FindBin qw ($Bin);
use lib "$Bin/lib";
use rawProbe;

use Devel::Size qw[ total_size ];
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
#use Parallel::ForkManager;


# get options
my $help = 0;
my $genomeFile;
my $targetFile;
my $kmerFile;
my $step = 5;
my $rep = 5;
my $amb = 1;
my $kmer = 17;
my $mode = "lowMem";
my $multiMer = 4;
my $homo = 0.7;
my $tmThr  = 37;
my $htmThr = 35;
my $dtmThr = 10;
my $thre=4;
my $uniqKmerCutoff=1;

GetOptions(
    'help' => \$help,
    'gs=s'   => \$genomeFile,   # genome sequence file
    't=s'    => \$targetFile,   # target sequence file
	'km=s'   => \$kmerFile,     # pre-calculated kmer frequency file
	'uk=i'   => \$uniqKmerCutoff, # uniq kmer cutoff, default is 1 for genomic, but if not use genomic for kmer profile, set accordingly
    's=i'    => \$step,         # step
    'r=i'    => \$rep,          # max length of repeat nucleotide (AAAA,TTTT...)
    'a=i'    => \$amb,          # max number of N or n
    'k=i'    => \$kmer,         # length of kmer
    'm=s'    => \$mode,         # running mode
    'n=i'    => \$multiMer,     # max number of multiple kmers in a probe
    'h=f'    => \$homo,         # max homologous in the genome
    'p=i'    => \$tmThr,        # min tm of probe
    'q=i'    => \$htmThr,       # max tm for hairpin (htm) of probe
    'd=i'    => \$dtmThr,       # min differences between tm and htm
	'e=i'    => \$thre          # number of threads
) or pod2usage(2);

pod2usage(1) if $help;

if(!(defined $genomeFile) || !(defined $targetFile)){
	print STDERR "please set a genome sequence(-gs) and a target sequence file (-t).\n";
	pod2usage(1);
}

if(!(-e $genomeFile)){
	print "\n***** $genomeFile doesn't exist. *****\n\n";
	pod2usage(1);
}


if(!(-e $targetFile)){
	print "\n***** $targetFile doesn't exist. *****\n\n";
	pod2usage(1);
}
	

print STDERR "genome sequencefile is $genomeFile\n
target sequence file is $targetFile\n
step is $step\n
max length of tendom repeat is $rep\n
max Nn in a probe is $amb\n
max homologous in the genome is $homo\n
kmer size is $kmer\n
running mode is $mode\n
running on $thre threads\n

thermodynamics are tm $tmThr, htm $htmThr, dtm is $dtmThr\n\n";


my $tt=time;
# generate candidate probe
print STDERR "generate candidate probe from genome in a step of $step bp \n";

my $obj;

$obj=rawProbe->new($targetFile,$genomeFile);
$obj->makeProbe($step,$rep,$amb);
print STDERR "make probe use ",time-$tt," seconds\n";

# generate kmer profile of candidate probe and whole genome

$obj->runKMC($kmer,$thre,$kmerFile,$uniqKmerCutoff);

# find kmer frequncy of each probe and get tm
print "analyze kmer frequency of each probe\n";
$mode="highMem";
$obj->kmerFreq($thre,$mode);



# check homologous

print STDERR "check homologous (max homologous is $homo)\n";

print STDERR "\tprepare reads\n";

$obj->filterKmer($multiMer);

print STDERR "\trun bwa\n";
my $bam=$obj->runBWA($thre);

print STDERR "\textract probes with less than $homo homologous\n";
my $homoFiltered=$obj->filterHomo($bam,$homo);

print STDERR "check tm of each probe\n";


my $filteredProbe=$obj->calTm($homoFiltered, $tmThr,$htmThr,$dtmThr,$thre);

print STDERR "\nall finished,use ",time-$tt," seconds\nfinal probe is stored in file $filteredProbe\n";

__END__

=head1 NAME

Chorus_perl - a perl impliment of Chorus, which allows design probes in low memory machine.

=head1 SYNOPSIS

perl Chorus_perl.pl -gs <genome file> -t <target sequence file>

Options:

    	-help              this help message
          -gs              genome sequence file, mandatory
           -t              target sequence file, mandatory
          -km              pre-calculated kmer profile
          -uk              uniq kmer cutoff, default is 1 for genomic, but if not use genomic for kmer profile, set accordingly
           -s              step of ajacent probe , default is every 5bp.
           -r              max length of repeat nucleotide (AAAA,TTTT...), default is 5
           -a              max number of N or n, default is 1
           -k              length of kmer, default is 17
           -m              running mode,canbe set to highMem and lowMem, default is highMem
           -n              max number of multiple kmers in a probe, default is 4
           -h              max homologous in the genome, default is 0.7
           -p              min tm of probe, default is 37
           -q              max tm for hairpin (htm) of probe, default is 35
           -d              min differences between tm and htm, default is 10
           -e              number of threads, default is 4


=head1 DESCRIPTION

Chorous_perl.pl is a perl impliment of Chorus (https://github.com/forrestzhang/Chorus). It do the samething as Chorous but introduce a low-memory running mode.

=cut


