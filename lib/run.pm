#
#===============================================================================
#
#         FILE: run.pm
#
#  DESCRIPTION:
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 03/17/2017 12:39:45 PM
#     REVISION: ---
#===============================================================================

package run;
use strict;
use warnings;

sub new {
    my $class = shift;
    my $self  = {};
    return bless $self, $class;
    return $self;
}    # ----------  end of subroutine new  ----------

sub runKMC {
    my $self = shift;

    my $kmer = shift;
    my $thre = shift;

	my $kmerProfile=shift;
	my $uniqKmerCutoff=shift;

	$uniqKmerCutoff=1 if !$uniqKmerCutoff;

    print STDERR "--------------calcuate kmer (k=$kmer)-----------------\n";

    my $targetFile = $self->{'targetFile'};
    my $genomeFile = $self->{'genome'};
    my $temProbe   = $self->{'rawProbe'};
    my $temDir     = $targetFile . ".temDir";

    my $canKmerProfile = $temProbe . ".k$kmer";
#    my $comKmer        = $kmerProfile . ".comm";
#    my $dump           = $comKmer . ".dump";

    if ( -d $temDir ) {
        print STDERR
            "\t******* directory $temDir already existed ***********\n";
    }
    else {
        print STDERR "******* make directory $temDir *******\n";
        mkdir "$temDir";
    }

    print STDERR "\t******* generate kmer profile *******\n";

	if(!$kmerProfile){
		$kmerProfile    = $genomeFile . ".k$kmer";
	    system "kmc -t$thre -k$kmer -ci1 -fm $genomeFile $kmerProfile $temDir";
	}
    system "kmc -t$thre -k$kmer -ci1 -fa $temProbe $canKmerProfile $temDir";

    print STDERR "\t******* get unique kmers *******\n";
	my $comKmer        = $kmerProfile . ".comm";
	my $dump           = $comKmer . ".dump";
    system
        "kmc_tools simple $canKmerProfile $kmerProfile -ci1 intersect $comKmer -ocright";
    system "kmc_tools transform $comKmer -ci1 -cx$uniqKmerCutoff dump -s $dump";


	print STDERR "\n\n\tuniq kmer out is $dump\n\n";

    $self->{'uniqKmer'} = $dump;
}

sub runBWA {
    my $self       = shift;
    my $threads    = shift;
    my $genomeFile = $self->{'genome'};
    my $temReads   = $self->{'temReads'} || die("no file for homology check\n");
    my $bamFile    = $temReads;
    my $filtered   = $temReads;
    $bamFile =~ s/reads$/bam/;
    $filtered =~ s/reads$/filtered/;
	print STDERR "run BWA\n\n";
	if(!(-e $genomeFile.".pac" || -e $genomeFile.".bwt" || -e $genomeFile.".ann" || -e $genomeFile.".amb" || -e $genomeFile.".sa")){
		print STDERR "\tbuilding bwa index\n";
		system "bwa index $genomeFile";
	}

	print STDERR "\trun alignment\n";
    system
        "bwa mem -O 0 -B 0 -E 0 -k 5 -t $threads $genomeFile $temReads | samtools view -bS - -o $bamFile";
    $self->{'homoBamFile'} = $bamFile;
	return $bamFile;
}    ## --- end sub runBWA
1;
