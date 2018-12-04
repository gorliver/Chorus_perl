# Chorus_perl
A perl implement of Chorus idea. Checkout https://github.com/forrestzhang/Chorus.
This project is just for learning purpose. Chorus_perl is by no means to take any credit. All the credits and citatation should refer to https://github.com/forrestzhang/Chorus.

# requirement

perl model:

    Inline::python
    Parallel::ForkManager
    BerkeleyDB

external tools:

    samtools (http://www.htslib.org/)
    bwa (https://sourceforge.net/projects/bio-bwa/files/)
    KMC3.0 (http://sun.aei.polsl.pl/kmc)


python:

    primer3-py (https://pypi.python.org/pypi/primer3-py)


# usage

    perl Chorus_perl.pl -gs <genome file> -t <target sequence file>

    Options:

            -help              this help message
              -gs              genome sequence file (while genome, including unanchored scaffolds/contigs), MANDATORY
               -t              target sequence file (could be one region or one chromosome. could be whole genome but not recommended),
                               MNADATORY
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
