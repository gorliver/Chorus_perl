package NewFa;

use 5.006;
use strict;
use warnings;

our $VERSION='0.01';


sub new {
	my $class=shift;
	my $fh=shift;
	my $self={FH=>$fh};
	bless $self, $class;
}

sub read_seq {
	my $self=shift;
	my $fh=$self->{FH};

	my $ll=0; # check if no new lines for parse	

	if(my $last=$self->{LAST}){ # header line 
		$last=~/^>(.*)/;
		$self->{DEF}=$1;
		
		my $seq;
		while(<$fh>){ # get the seq
			$ll++;
			last if /^>/;
			chomp;
			$seq.=$_;
		}
		$self->{LAST}=$_;
		$self->{SEQ}=$seq;
		$self->{LEN}=length $seq;
		return 0 if $ll==0;
		return $self;
	}
	
	
	
	while(<$fh>){
		$ll++;
		if(/^>(.*)/){
			my $seq;
			$self->{DEF}=$1; # get header.
			while(<$fh>){
				last if /^>/;
				chomp;
				$seq.=$_;
			}
			$self->{LAST}=$_; # record the next header (or empty if reach the end of the file)
			$seq=~s/\s+//g;
			$self->{SEQ}=$seq; # get seq
			$self->{LEN}=length $seq;
			last;
		}
	}
	return 0 if $ll==0;
	return $self;
}

sub def {
	my $self=shift;
	$self->{DEF};
}
sub seq {
	my $self=shift;
	$self->{SEQ};
}

sub len {
	my $self=shift;
	$self->{LEN};
}

1;

__END__

=head1 NAME

 NewFa;

=head1 SYNOPSIS

 use NewFa;


 open $ff,$ARGV[0] or die;
 my $fasta=newFa->new($ff);
 while(my $one=$fasta->read_seq){
 	print $one->def,"\n";
 	print $one->seq,"\n";
 }


=head1 DESCRIPTION

newFa is a package to deal with stupid multi-fasta files, extract the header (/>(\S+)) and the seq pairs.

