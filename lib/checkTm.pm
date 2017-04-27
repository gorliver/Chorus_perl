
#
#===============================================================================
#
#         FILE: tmCheck.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 03/20/2017 10:56:46 AM
#     REVISION: ---
#===============================================================================
package checkTm;
use strict;
use warnings;
use Parallel::ForkManager;
use parent 'kmerFilter';

sub new {
	my $class=shift;
	my $self={};
	return bless $self,$class;
}

sub calTm {
	my $self=shift;
	my $probe=shift;
	my $tmThr=shift;
	my $htmThr=shift;
	my $dtmThr=shift;
	my $threads=shift;

	my $filtered=$probe.".filtered";

	print STDERR "\tanalyzing $probe in $threads threads\n";
	my @temOutFiles=@{$self->splitFile($probe,$threads)};

	my @temOutTM;
	my $tmpm=Parallel::ForkManager->new($threads);
	$tmpm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
				$data_structure_reference )
				= @_;

			if ( defined($data_structure_reference) ) {
				push @temOutTM, $$data_structure_reference;
			}
			else {
				print STDERR "process $pid did not return a result, it probably failed.\n\n";
				exit;
			}
		}
	);

	TM:
	foreach my $oneTemFile(@temOutFiles){
		my $pid=$tmpm->start and next TM;
		
		my $outTem=$oneTemFile . ".out";
		open FF,$oneTemFile or die;
		open OO,"+>$outTem" or die;
		while(<FF>){
			my @tem=split;
			my ($tm,$htm)=primer3_filter($tem[9]);
			my @arr=split /\_/,$tem[0];
			if($tm>$tmThr && $htm<$htmThr && $tm-$htm>$dtmThr){
				print OO "$tem[9]\t$arr[0]\t$arr[1]\t",$arr[1]+44,"\t$arr[2]\t$arr[3]\t\t$tm\t$htm\n";
			}
		}
		close FF;
		close OO;

		$tmpm->finish(0,\$outTem);
	}

	print STDERR "\tall tm job submitted\n";
	$tmpm->wait_all_children;
	print STDERR "\tall tm job done, got @temOutTM, merge results\n";
	system "cat @temOutTM | sort -k2,2 -k3,3n > $filtered";
	print STDERR "\tclean tm jobs\n";
	system "rm @temOutTM";
	system "rm @temOutFiles";
	system "rm -rf _Inline/";

	my ($tp)=split /\s+/,`wc -l $filtered`;
	my ($c2)=split /\s+/,`wc -l $probe`;

	print STDERR "\tgot $c2 probes, of which $tp probe satisfied the thermodynamics (min Tm is $tmThr, max hTm is $htmThr, min(Tm-hTm) is $dtmThr)\n";

	return $filtered;
}


use Inline Python => <<'END_PYTHON';

import primer3

def primer3_filter(sequence, mintm=37, maxhtm=35, dtm=10):

   tm = primer3.calcTm(sequence)
   htm = primer3.calcHairpinTm(sequence)

   return (tm,htm)


END_PYTHON


1;

