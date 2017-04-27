#
#===============================================================================
#
#         FILE: rawProbe.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 03/17/2017 12:12:20 PM
#     REVISION: ---
#===============================================================================
package rawProbe;

use strict;
use warnings;
use parent 'run';
use parent 'kmerFilter';
use parent 'checkTm';
use NewFa;

sub new{
	my $class=shift;
	my $targetFile=shift;
	my $genomeFile=shift;
	my $rawProbe=$targetFile.".tem.fa";

	my $self={
		'targetFile'=>$targetFile,
		'genome'=>$genomeFile,
		'rawProbe'=>$rawProbe
	};
	return bless $self, $class;
}

sub makeProbe{
	my $self=shift;

	my $step=shift;
	my $rep=shift;
	my $amb=shift;

	my $targetFile=$self->{'targetFile'};
	my $temProbe=$self->{'rawProbe'};

	my $g1=0;
	my $c1=0;
	open my $tf,$targetFile or die ("cannot find $targetFile\n");
	open O,"+>$temProbe" or die ("cannot create $temProbe\n");
	my $fasta=NewFa->new($tf);
	while(my $one=$fasta->read_seq){
		my $name=$one->def;
		my $seq=$one->seq;
		my $len=$one->len;
		print STDERR "length is $len\n";
		my $st=0;
		my $repeatLen;
		my $repeatEd;
		my $repeatSt;
		my $nonRepLen;
		while($seq=~/(A{$rep,}|T{$rep,}|C{$rep,}|G{$rep,})/g ) {
			$repeatLen = length $1;
			$repeatEd  = pos($seq);
			$repeatSt  = $repeatEd - $repeatLen;
			$nonRepLen = $repeatSt - $st;
			$g1+=$repeatLen;
			$g1+=$nonRepLen;
			my $nr = substr( $seq, $st, $repeatSt - $st );
			for ( my $x = 0; $x <= $nonRepLen - 45; $x += 5 ) {
				$c1++;
				my $subSeq=substr( $seq, $st + $x, 45 );
				next if $subSeq=~/[Nn]{$amb,}/;
				print O ">$name\_", $st + $x, "\n", substr( $seq, $st + $x, 45 ), "\n";
				print STDERR "generate $c1 for ",$g1+1," bp\n" if !($c1%1000000);
			}
			$st = $repeatEd;
		}
		$nonRepLen=$repeatEd-$st;
		for ( my $x = 0; $x <= $nonRepLen - 45; $x += 5 ) {
			print O ">$name\_", $st + $x, "\n", substr( $seq, $st + $x, 45 ), "\n";
		}
	}
	close $tf;
	close O;
	print STDERR "\tanalyzed $g1 bp, got $c1 candidate probes\n";
}


sub skip {
	my $self=shift;

	my $ss=shift;

	my @arr=split / /,$ss;
	my $nn=@arr;
	if($nn<=1){
		return 45;
	}
	else{
		return 45-(length $arr[-1]);
	}
}

1;
