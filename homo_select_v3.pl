#!/usr/bin/perl -w

my %target;
my %notCore;
my $lastChr=0;
my $lastPos=0;
open F, "samtools view $ARGV[0] | " or die;
while(<F>){
	my @tem=split;

	/AS:i:(\d+)/;
	next if $1<40;
	my $score=0;

	my $mq=$tem[4];
	my ($mis)=/NM:i:(\d+)/;

	my @indel=$tem[5]=~/(\d+)[ID]/g;

	my $indelSum=0;
	grep {$indelSum+=$_} @indel;

	my @arr=split /\_/,$tem[0];

	next if $indelSum>0; # no mismatch

	next if $mq<50; # unique in target genome

	if($arr[0] ne $lastChr){
		$lastChr=$arr[0];
		$lastPos=$arr[1];
		redo;
	}
	if($mis==0){ # perfect
		if($arr[1]-$lastPos>5){
			$target{$arr[0]}{$arr[1]}=[$tem[9],$mq,$mis];
			$lastPos=$arr[1]+44;
		}
	}
	else{# less perfect
		if($arr[1]-$lastPos>5){
			$target{$arr[0]}{$arr[1]}=[$tem[9],$mq,$mis];
			$lastPos=$arr[1]+44;
		}
	}

	$lastChr=$arr[0];


}
close F;

print STDERR "reading bam done\n";

open G,$ARGV[1] or die;
while(<G>){
	my @tem=split;
	$notCore{$tem[0]}{$tem[1]}=[$tem[3],100,100] if !$notCore{$tem[0]}{$tem[1]};
}
close G;

print STDERR "reading file done\n";
my @chr=keys %notCore;
my @chr2=keys %target;

print STDERR "chrs are @chr, @chr2\n";
foreach my $nChr(@chr){
	next if !$target{$nChr};
	print STDERR "analysis $nChr\n";

	my @notCorePos=sort {$a<=>$b} keys %{$notCore{$nChr}};
	my @work=@notCorePos;
	my @corePos=sort {$a<=>$b} keys %{$target{$nChr}};

	my $nnc=@corePos;
	my $nnnc=@notCorePos;
	print STDERR  "got $nnc core and $nnnc notcore\n";
	my $numCore=@corePos;

	my $restPrb=27392-$numCore;

	print STDERR "have $numCore perfect probe, need assign $restPrb to poor regions\n";
	# get distribution of core probe each 100kb window
	
	my %distr;
	foreach my $onePrb(@corePos){
		my $id=int($onePrb/100000);
		$distr{$id}++;
	}

	# assign not_ core probe to the poor windows

	my @zeroWin;
	my $poorSt=0;
	for(my $x=0;$x<int($notCorePos[-1]/100000);$x++){
		if(!$distr{$x}){
			push @zeroWin,$x;
		}
	}

	my $zeroWinNum=@zeroWin;
	print STDERR " got $zeroWinNum window without any probe\n";

	my $assignNum=int($restPrb/$zeroWinNum+0.5);
	$assignNum=1 if $assignNum==0;
	print STDERR " assign probe to zero-window, each window will have $assignNum\n";

	my @posSortByMisThenOrder=sort {$notCore{$nChr}{$a}->[2] <=> $notCore{$nChr}{$b}->[2] or
						$a<=>$b
						} keys %{$notCore{$nChr}};


	my $ncNum=@posSortByMisThenOrder;
	print STDERR "we have $ncNum to add \n";
	my $total=0;
	# set mach to avoid over represent at poor region
	my $max=int(27392/($notCorePos[-1]/100000));
	$max=int($max*1.2);
#	$max=$assignNum;
	($total,$restPrb)=assignProbe(\@notCorePos,\@corePos,\@zeroWin,$nChr,$restPrb,$max,$total);
	print STDERR "assigned ",$total-1," probe\n";
	
	if($restPrb>0){
		my @nonZeroWin=sort {$a<=>$b} keys %distr;
		my $assignNum2=int($restPrb/@nonZeroWin+0.5);
		$assignNum2=1 if $assignNum2==0;
		print STDERR "not enough, assign rest $restPrb randomly, each will have $assignNum2\n";
		($total,$restPrb)=assignProbe(\@notCorePos,\@corePos,\@nonZeroWin,$nChr,$restPrb,$assignNum2,$total);

		print STDERR "second assigned ",$total-1," probe\n";
	}
}

print STDERR "assign done\n";

foreach my $nc(keys %target){
	my @arr=sort {$a<=>$b} keys %{$target{$nc}};

	foreach my $np(@arr){
		print "$nc\t$np\t",$np+44,"\t",join "\t",@{$target{$nc}{$np}},"\n";
	}
}


sub getBest{# inserted candidate probe to target probe set
	my $candiMisRef=shift;
	my $chr=shift;
	my $asNum=shift;
	my $good = shift;
	my $maxPrb=shift;
	my $res=shift;
	my $added=0;
	my $idx=0;
#	print "got $maxPrb, $res\n";
	my $switch = $maxPrb;#>$asNum ? $maxPrb : $asNum;
	while(1){
		if(!@{$good}){
			push @{$good},$$candiMisRef[$idx];
			$target{$chr}{$$candiMisRef[$idx]}=$notCore{$chr}{$$candiMisRef[$idx]};
			$added++;
			$res--;
		}
		else{
			for(my $i=0;$i<@{$good};$i++){
				if($$good[$i]-($$candiMisRef[$idx]+44)>5 ){ #pre
					if($i==0 || $$candiMisRef[$idx]-($$good[$i-1]+44)>5){
						splice @{$good},$i,0,$$candiMisRef[$idx];
						$added++;
						$res--;
						$target{$chr}{$$candiMisRef[$idx]}=$notCore{$chr}{$$candiMisRef[$idx]};
						last;
					}
				}
				elsif($$candiMisRef[$idx]-($$good[$i]+44)>5){#suc
					if(!$$good[$i+1] || $$good[$i+1]-($$candiMisRef[$idx]+44)>5){
						splice @{$good},$i+1,0,$$candiMisRef[$idx];
						$added++;
						$res--;
						$target{$chr}{$$candiMisRef[$idx]}=$notCore{$chr}{$$candiMisRef[$idx]};
						last;
					}
				}
			}

			$idx++;
			last if !$$candiMisRef[$idx];
#			if($added>=$asNum){
			if($added>=$switch || $res<=0){
				last;
			}
		}
	}

	return ($good,$added,$res);
}

sub assignProbe{
	my $probeRef=shift;
	my $targetRef=shift;
	my $windowRef=shift;
	my $chr=shift;
	my $restNum=shift;
	my $maxNum=shift;
	my $total=shift;
#	print "this max is $maxNum\n";
	my @work=@{$probeRef};
	my @tar=@{$targetRef};

	my $winNum=@{$windowRef};

	my $asNum=int($restNum/$winNum+0.5);

	my $x=0;
    while($winNum){
        my $st=$$windowRef[$x]*100000;
        my $ed=$st+99999;

        my @candi;
		my @tarReg;
      	while($work[0]<$ed){
			push @candi,$work[0] if $work[0]>$st;
			shift @work;
		   	last if !$work[0];
		}

		while($tar[0]<$ed){
			push @tarReg,$tar[0] if $tar[0]>$st;
			shift @tar;
			last if !$tar[0];
		}

        my $tt=@candi;

        my @candiMis=sort {$notCore{$chr}{$a}->[2] <=> $notCore{$chr}{$b}->[2]} @candi;

        my @good;
		my $goodRef;
		my $added;
        ($goodRef,$added,$restNum)=getBest(\@candiMis,$chr,$asNum,\@tarReg,$maxNum,$restNum);

		$total+=$added;

        $winNum--;
        $x++;
        last if $restNum<=0 || $winNum==0 || !$$windowRef[$x];
#        $asNum=int($restNum/$winNum+0.5);

        @good=@{$goodRef};
        my $ttt=@good;
   }
   return $total,$restNum;
}

