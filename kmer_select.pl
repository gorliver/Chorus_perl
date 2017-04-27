#!/usr/bin/perl -w

my $str="+";
#my $lastf;
#my $lastr;
open F,$ARGV[0] or die;
while(<F>){
	exit if !$_;
	my @tem=split;
	next if $tem[-1] eq "false";

#	print "it is $_\n";
	if($str eq "-"){
		$tem[3]=reverse $tem[3];
		$tem[3]=~tr/ATCG/TAGC/;
#		$lastr=$tem[2];
	}
#	else{
#		$lastf=$tem[2];
#	}

	print join "\t",(@tem,$str);
	print "\n";

	exit if !$_;
	while(<F>){
		exit if !$_;
		my @new=split;
#		print "\tcheck $_\n";
		next if $new[-1] eq "false";
		
		if($new[1]<$tem[2]+5 && $new[1]>$tem[1]+15){
#			print "\t\toverlap\n";
			$str= $str eq "+" ? "-" : "+";
			last;
		}
		if($new[1]>=$tem[2]+5){
#			print "\t\tnonoverlap\n";
			$str= $str eq "+" ? "-" : "+";
			last;
		}
	}

	redo;
}
close F;
