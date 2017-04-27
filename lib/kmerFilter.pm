#
#===============================================================================
#
#         FILE: kmerFilter.pm
#
#  DESCRIPTION:
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 03/17/2017 02:04:40 PM
#     REVISION: ---
#===============================================================================
package kmerFilter;
use strict;
use warnings;
use Parallel::ForkManager;
use BerkeleyDB;
use Devel::Size qw[ total_size ];

sub new {
    my $class = shift;
    my $self  = {};
    return bless $self, $class;
}

sub filterHomo {
	my $self=shift;
	my $bamFile=shift;
	my $homo=shift;
	my $homoNum=int(45*$homo);
	my $temOut=$bamFile;
	$temOut=~s/bam/homo/;
	open B,"samtools view $bamFile | " or die;
	open TO,"+>$temOut" or die;
	while(<B>){
		next if /^\@/;
		my @tem=split;
		/XS:i:(\d+)/;
		if($1<=$homoNum){
			print TO $_;
		}
	}
	close B;
	close TO;

	return $temOut;

}
sub splitFile {
    my $self    = shift;
    my $file    = shift;
    my $threads = shift;
    my @files;

    my ($lineNum) = split /\s+/, `wc -l $file`;
    print STDERR "\t\t***** have $lineNum lines\n";
    my $oneNum = int( $lineNum / $threads + 1 );
    $oneNum++ if $oneNum % 2;
    print STDERR "\t\t***** have $oneNum lines for each file\n";

    my $id  = 0;
    my $out = $file . ".$id";
    push @files, $out;
    open TF, $file    or die;
    open O,  "+>$out" or die;
    my $nn = 0;

    while (<TF>) {
        $nn++;
        if ( $nn <= $oneNum ) {
            print O $_;
        }
        else {
            close O;
            $id++;
            $out = $file . ".$id";
            $nn  = 0;
            push @files, $out;
            open O, "+>$out" or die;
            redo if $_;
        }
    }
    close TF;
    close O;

    print STDERR "\t\t$file were devided into $threads files\n";

    return \@files;
}
sub filterKmer {
    my $self     = shift;
    my $multiMer = shift;
    my $kmerOut  = $self->{'kmerCount'} || die("no result from kmerFreq\n");
    my $temReads = $kmerOut . ".m$multiMer\.reads";
    my $c1       = 0;
    open K, $kmerOut
        or die("cannot open kmerFreq file. kmer analysis probably failed.\n");
    open TR, "+>$temReads" or die("cannot create $temReads\n");

    while (<K>) {
        my @tem = split;
        if ( $tem[4] <= $multiMer ) {
            print TR ">$tem[0]_$tem[3]_$tem[4]\n$tem[1]\n";
            $c1++;
        }
    }
    close K;
    close TR;

    print STDERR
        "\tgot $c1 probes contained less than $multiMer kmers which are appeared multiple times in the genome\n";
	$self->{'temReads'}=$temReads;
}

sub kmerFreq {
    my $self    = shift;
    my $threads = shift;
    my $mode    = shift;
    print STDERR
        "-------------- analyze kmer frequency of each probe ----------------\n";

    if ( $mode eq "highMem" ) {
        print STDERR "\trunning in highMem mode\n";
    }
    elsif ( $mode eq "lowMem" ) {
        print STDERR "\trunning in lowMem mode\n";
    }
    else {
        print STDERR "\tmode has to be 'highMem' or 'lowMem'\n";
        exit;
    }

    my $dumpFile = $self->{'uniqKmer'}
        || die("cannot find uniqKmer. kmer analysis probably failed.\n");
    my $temProbeFile = $self->{'rawProbe'}
        || die("cannot find uniqKmer. kmer analysis probably fai    led.\n");

    my $output = $temProbeFile . ".filtered";

    # split temProbeFile
	my @files=@{$self->splitFile($temProbeFile,$threads)};
	my @filesOut;

	print STDERR "\t\t$temProbeFile were devided into $threads files\n";

    # make database
    if ( $mode eq "highMem" ) {
        print STDERR "\tmake database\n";
        my $startPoint = time;
        my %hash;
        my $aa = 0;
        open DUMP, $dumpFile or die;
        while (<DUMP>) {
            $aa++;
            print STDERR "done $aa\n" if !( $aa % 5000000 );
            my @tem = split;
            $hash{ $tem[0] } = $tem[1];
        }
        close DUMP;

        my $tt = time - $startPoint;
        printf STDERR "running time: %d seconds\n", $tt;
        printf STDERR "Memory: %f MB\n", total_size( \%hash ) / 1024**2;

        print STDERR
            "\t\twe got @files, do the analysis in $threads threads\n";
        my $kmerpm = Parallel::ForkManager->new($threads);

        $kmerpm->run_on_finish(
            sub {
                my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                    $data_structure_reference )
                    = @_;
                if ( defined($data_structure_reference) ) {
                    push @filesOut, $$data_structure_reference;
                }
                else {
                    print STDERR
                        "process $pid did not return a result, it probably failed.\n\n";
                    exit;
                }
            }
        );

    KMER:
        foreach my $oneFile (@files) {
            my $pid = $kmerpm->start and next KMER;

            my $outTm = $oneFile . ".out";
            my $cc    = 0;
            open OF, $oneFile   or die;
            open O,  "+>$outTm" or die;

            my $name;
            my $seq;
            while (<OF>) {
                if (/^>(\S+)/) {
                    $name = $1;
                    next;
                }
                chomp;
                $seq = $_;

                $cc++;
                $tt = time - $startPoint;
#                print STDERR "$tt done $cc\n" if !( $aa % 1000000 );
                my $len        = length $seq;
                my $kcount     = 0;
                my $kcountm    = 0;
                my $kcountmid  = 0;
                my $kcountmidm = 0;
                my $keep       = "true";
                my $uc         = 0;
                my $mc         = 0;
                my $sum        = 0;
                for ( my $x = 0; $x <= $len - 17; $x++ ) {
                    my $sub = substr( $seq, $x, 17 );
                    my $rev = reverse $sub;
                    $rev =~ tr/ATCG/TAGC/;
                    if ( $hash{$sub} || $hash{$rev} ) {
                        $uc++;
                    }
                    else {
                        $mc++;
                    }
                }
                my $freq = $uc + $mc;
                $keep = "false" if $freq < 5 || $freq > 70;

                print O "$name\t$seq\t$keep\t$uc\t$mc\n";    #\t$tm\t$htm\n";
            }

            close OF;
            close O;

            $kmerpm->finish( 0, \$outTm );
        }

        print STDERR "\tall kmer job submitted\n";
        $kmerpm->wait_all_children;
        print STDERR "\tall kmer job done, merge the results\n";
        system "cat @filesOut > $output";

        print STDERR "\tclean kmer job\n";
        system "rm @files";
        system "rm @filesOut";
        print STDERR "\tuse ", time - $startPoint, " seconds\n";

        undef %hash;
    }
    else {
        print STDERR "\tmake database\n";
        my $startPoint = time;
        my $filename   = "database";
        unlink $filename;
        my $dbh = new BerkeleyDB::Hash
            -Filename => $filename,
            -Flags    => DB_CREATE
            or die "Cannot open file $filename : $! $BerkeleyDB::Error\n";

        my $aa = 0;
        open DUMP, $dumpFile or die;
        while (<DUMP>) {
            $aa++;
            print STDERR "done $aa\n" if !( $aa % 5000000 );
            my @tem = split;
            $dbh->db_put( $tem[0], $tem[1] );
        }
        close DUMP;

        my $tt = time - $startPoint;
        printf STDERR "running time: %d seconds\n", $tt;
        printf STDERR "Memory: %f MB\n", total_size($dbh) / 1024**2;

        my @filesOut;
        print STDERR
            "\t\twe got @files, do the analysis in $threads threads\n";
        my $kmerpm = Parallel::ForkManager->new($threads);

        $kmerpm->run_on_finish(
            sub {
                my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
                    $data_structure_reference )
                    = @_;
                if ( defined($data_structure_reference) ) {
                    push @filesOut, $$data_structure_reference;
                }
                else {
                    print STDERR
                        "process $pid did not return a result, it probably failed.\n\n";
                    exit;
                }
            }
        );

    KMER:
        foreach my $oneFile (@files) {
            my $pid   = $kmerpm->start and next KMER;
            my $outTm = $oneFile . ".out";
            my $cc    = 0;
            open OF, $oneFile   or die;
            open O,  "+>$outTm" or die;
            my $name;
            my $seq;
            while (<OF>) {

                if (/^>(\S+)/) {
                    $name = $1;
                    next;
                }

                chomp;
                $seq = $_;

                $cc++;
                if ( !( $aa % 1000000 ) ) {
                    $tt = time - $startPoint;
                    print STDERR "$tt done $cc\n";
                }
                my $len        = length $seq;
                my $kcount     = 0;
                my $kcountm    = 0;
                my $kcountmid  = 0;
                my $kcountmidm = 0;
                my $keep       = "true";
                my $uc         = 0;
                my $mc         = 0;
                my $sum        = 0;
                for ( my $x = 0; $x <= $len - 17; $x++ ) {
                    my $sub = substr( $seq, $x, 17 );
                    my $rev = reverse $sub;
                    $rev =~ tr/ATCG/TAGC/;
                    if (   !$dbh->db_exists($sub)
                        || !$dbh->db_exists($rev) )
                    {    # if not exit will return DB_NOTFOUND
                        $uc++;
                    }
                    else {
                        $mc++;
                    }
                }
                my $freq = $uc + $mc;
                $keep = "false" if $freq < 5 || $freq > 70;

                print O "$name\t$seq\t$keep\t$uc\t$mc\n";    #\t$tm\t$htm\n";
            }

            close OF;
            close O;

            $kmerpm->finish( 0, \$outTm );
        }

        print STDERR "\tall kmer job submitted\n";
        $kmerpm->wait_all_children;
        print STDERR "\tall kmer job done, merge the results\n";
        system "cat @filesOut > $output";

        print STDERR "\tclean kmer job\n";
        system "rm @files";
        system "rm @filesOut";
        print STDERR "\tuse ", time - $startPoint, " seconds\n";

        undef $dbh;
    }
    $self->{'kmerCount'} = $output;

}    ## --- end sub kmerFreq

1;
