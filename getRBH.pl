#!/usr/bin/perl
## a little program to extract orthologs from a homologs file
## using the reciprocal best hits definition
## an example line is printed when no arguments are given
use strict;
use Getopt::Long;
use Pod::Text;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);
my $columns = qx(tput cols);
chomp($columns);
my $width = $columns >= 100 ? 80 : $columns - 3;
my $parser = Pod::Text->new (sentence => 0, width => $width, margin => 1);

my $ownName = $0;
$ownName =~ s{.*/}{};

my %url = (
    "blastp"  => 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/',
    "diamond" => 'https://github.com/bbuchfink/diamond/releases',
    "mmseqs"  => 'https://github.com/soedinglab/MMseqs2/releases',
    "lastal"  => 'http://last.cbrc.jp/',
);
my @pwProgs = sort keys %url;
my $matchProg = join("|",@pwProgs);
my ($refavailable,$refmissing,$missMsg) = checkSoftware();

###### defaults, etc
my $cov1       = 50;
my $cov2       = 99;
my $defCov     = 60;
my $defAln     = "F";
my $defPW      = 'blastp';
###### assignable:
my @queries    = ();
my @against    = ();
my $faaDir     = '';
my $rbhDir     = "RBH";
my $minCov     = $defCov;
my $maxOverlap = 0.1; # maximum overlap in fused proteins/genes
my $maxEvalue  = 1e-6;
my $alnSeqs    = $defAln;
my $pwProg     = $defPW;
my $pwDir      = "compRuns";
my $sensitive  = 'S';
my $runRBH     = "T";
my $keepold    = 'T';
my $cpus       = 1;
my $topMatch   = 0;
my $minTop     = 50; # 30 was enough for lastsl in one example
my $matchRatio = 7;

### check if there's more than one processor or assume there's 1.
my $cpu_count
    = qx(sysctl -a | grep 'cpu.thread_count')
    =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
    =~ m{\.max-threads\s+=\s+(\d+)} ? $1
    : 1;

my $podUsage
    = qq(=pod\n\n)
    . qq(=head1 NAME\n\n)
    . qq($ownName - compare proteomes and extract reciprocal best hits\n\n)
    . qq(=head1 SYNOPSIS\n\n)
    . qq($ownName -q [filename1.faa] -t [filename2.faa] [options]\n\n)
    . qq(=head1 EXAMPLES\n\n)
    . qq($ownName -q Genomes/GCF_000005845.faa.gz -t Genomes/GCF_000009045.faa.gz\n\n)
    . qq($ownName -d Genomes -o AllvsAllRBH\n\n)
    . qq(=head1 OPTIONS\n\n)
    . qq(=over\n\n)
    . qq(=item B<-q>\n\n)
    . qq(query fasta file(s) or directory with fasta files, required.\n)
    . qq(files can be compressed with gzip or bzip2.\n\n)
    . qq(=item B<-t>\n\n)
    . qq(target fasta file(s) or directory with fasta files, required.\n)
    . qq(files can be compressed with gzip or bzip2.\n\n)
    . qq(=item B<-d>\n\n)
    . qq(directory with fasta files for all-vs-all comparisons, no default.\n)
    . qq(If set -q and -t will be ignored.\n\n)
    . qq(=item B<-o>\n\n)
    . qq(directory for reciprocal best hits, default: RBH.\n\n)
    . qq(=item B<-p>\n\n)
    . qq(program for pairwise comparisons [$matchProg], default: $defPW\n\n)
    . qq(=item B<-s>\n\n)
    . qq(sensitivity (only works for diamond), default '$sensitive':\n\n)
    . qq(=over\n\n)
    . qq(=item B<F:>\n\nfast, low sensitivity\n\n)
    . qq(=item B<S:>\n\nsensitive mode\n\n)
    . qq(=item B<X:>\n\nmore sensitive mode\n\n)
    . qq(=back\n\n)
    . qq(=item B<-m>\n\n)
    . qq(directory for $matchProg results, default compRuns.\n\n)
    . qq(=item B<-c>\n\n)
    . qq(minimum coverage of shortest sequence [$cov1 - $cov2], default $defCov.\n\n)
    . qq(=item B<-f>\n\n)
    . qq(maximum overlap between fused sequences [0 - 0.5], default 0.1\n)
    . qq((10%) of shortest sequence.\n\n)
    . qq(=item B<-e>\n\n)
    . qq(maximum e-value, default $maxEvalue\n\n)
    . qq(=item B<-a>\n\n)
    . qq(include aligned sequences in comparison results [T/F], default $defAln. (not available in lastal)\n\n)
    . qq(=item B<-r>\n\n)
    . qq(find RBH [T|F], default 'T', if 'F' the program will only run the\n)
    . qq($matchProg comparisons\n\n)
    . qq(=item B<-k>\n\n)
    . qq(keep prior results [T|R|F], default 'T':\n\n)
    . qq(=over\n\n)
    . qq(=item B<F:>\n\nrepeat both alignments and gathering of RBH\n\n)
    . qq(=item B<R:>\n\nonly repeat gathering of RBHs\n\n)
    . qq(=item B<T:>\n\nkeep prior alignments and prior RBHs\n\n)
    . qq(=back\n\n)
    . qq(=item B<-z>\n\nmaximum number of targets to find with $matchProg\n)
    . qq(minimum of $minTop. Defaults to 1/${matchRatio}th of the\n)
    . qq(sequences in the target file\n\n)
    . qq(=item B<-x>\n\nnumber of CPUs to use, default: 1 (max: $cpu_count)\n\n)
    . qq(=back\n\n)
    . qq(=head1 REQUIREMENTS\n\n)
    . qq(This program requires either appropriately formatted\n)
    . qq(comparison results for the query/target genomes, or\n)
    . qq($matchProg to compare protein sequences\n\n)
    . qq(=cut\n\n)
    ;

if( length($missMsg) > 0 ) {
    $podUsage
        .= qq(=head1 WARNING\n\n)
        .  qq(The following software is missing:\n\n)
        .  $missMsg;
}

GetOptions(
    "q=s{,}" => \@queries,
    "t=s{,}" => \@against,
    "d=s"    => \$faaDir,
    "p=s"    => \$pwProg,
    "m=s"    => \$pwDir,
    "s=s"    => \$sensitive,
    "o=s"    => \$rbhDir,
    "c=s"    => \$minCov,
    "f=s"    => \$maxOverlap,
    "e=s"    => \$maxEvalue,
    "a=s"    => \$alnSeqs,
    "r=s"    => \$runRBH,
    "k=s"    => \$keepold,
    "z=s"    => \$topMatch,
    "x=i"    => \$cpus,
) or podhelp();

### check if the selected program works
if( $refmissing->{"$pwProg"} ) {
    podhelp("$pwProg is missing");
}

### check if working with all-vs-all directory
my $allvsall = length($faaDir) > 0 ? 1 : 0;
if( $allvsall == 0 ) {
    if( !$queries[0] || !$against[0] ) {
        podhelp("if not -d I need both query and target files");
    }
    else {
        ####### check queries
        my %qseen = ();
        @queries = grep { not $qseen{"$_"}++ } @queries;
        my $qtoFind = @queries;
        if( -d $queries[0] && $qtoFind == 1 ) {
            @queries = findFaaFiles("$queries[0]","query");
        }
        else {
            @queries = findWorkFiles("query",@queries);
        }
        ###### check targets
        my %tseen = ();
        @against = grep { not $tseen{"$_"}++ } @against;
        my $stoFind = @against;
        if( -d $against[0] && $stoFind == 1 ) {
            @against = findFaaFiles("$against[0]","target");
        }
        else {
            @against = findWorkFiles("target",@against);
        }
    }
}
else {
    if ( -d "$faaDir" ) {
        if( @queries = findFaaFiles("$faaDir","all vs all") ) {
            @against = @queries;
        }
    }
    else {
        die "$faaDir doesn't exist or is not a directory\n\n";
    }
}

my $pwProg = $pwProg =~ m{^($matchProg)$} ? lc($1) : $defPW;
if( $pwProg eq "lastal" && $alnSeqs eq "T" ) {
    my $help
        = qq($pwProg does not save sequence alignments, if you need\n)
        . qq(the aligned sequences choose another program for pairwise\n)
        . qq(comparisons, otherwise don't use '-a T'\n\n);
    die "$help";
}
print "will compare proteins with $pwProg\n";
my $sensitive = $sensitive =~ m{^(F|S|X)$}i ? uc($1) : 'S';
my $dmdmode = '';
if( $pwProg eq "diamond" ) {
    $dmdmode
        = $sensitive eq 'F' ? "fast"
        : $sensitive eq 'S' ? "sensitive"
        : $sensitive eq 'X' ? "more-sensitive"
        : 'na';
    print "will run diamond in '$dmdmode' mode\n";
}
my $minCov = $minCov >= $cov1 && $minCov <= $cov2 ? $minCov : $defCov;
print "minimum coverage of shortest sequence: $minCov\n";
my $maxOverlap
    = $maxOverlap >= 0 && $maxOverlap <= 0.5
    ? $maxOverlap : 0.1;
print "maximum overlap for fusions: $maxOverlap\n";
my $alnSeqs = $alnSeqs =~ m{^(T|F)$}i ? uc($1) : $defAln;
print "include aligned seqs in blast results: $alnSeqs\n";
my $cpus = $cpus > 0 && $cpus <= $cpu_count ? $cpus : 1;
print "uing $cpus cpu threads\n";
my $runRBH = $runRBH =~ m{^(T|F)$}i ? uc($1) : 'T';
if( $runRBH eq 'F' ) {
    print "only pairwise comparisons will be run\n";
}
my $keepold = $keepold =~ m{^(T|F|R)$}i ? uc($1) : 'T';
if( $keepold eq 'T' ) {
    print "will keep all previous results\n";
}
else {
    if( $keepold eq 'R' ) {
        print "will gather RBHs even if a RBH file exists\n";
    }
    else {
        print "will make sequence comparisons and gather RBHs\n"
            . " even if prior results exists\n";
    }
}
#### directories where to find files:
#### temporary working directory:
my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");

##### headingd for pairwise comparisons:
my @pwHeading = qw(
                      Query
                      Target
                      eValue
                      bitScore
                      pident
                      qStart
                      qEnd
                      qLength
                      tStart
                      tEnd
                      tLength
              );
if( $alnSeqs eq "T" ) {
    push(@pwHeading,"qSeq","tSeq");
}

my @rbhHeading = qw(
                       Query
                       Target
                       eValue
                       bitScore
                       qStart
                       qEnd
                       qCover
                       tStart
                       tEnd
                       tCover
                       Qualifier
               );

##### the table format below is more informative than the default
##### works for both blastp and diamond
my @pwTbl = qw(
                  qseqid
                  sseqid
                  evalue
                  bitscore
                  pident
                  qstart
                  qend
                  qlen
                  sstart
                  send
                  slen
          );
##### add alignments to blast results table?
if( $alnSeqs eq "T" ) {
    push(@pwTbl,"qseq","sseq");
}
my $pwTbl = join(" ","6",@pwTbl);

##### table format for mmseqs:
my @mmfields = qw(
                     query
                     target
                     evalue
                     bits
                     pident
                     qstart
                     qend
                     qlen
                     tstart
                     tend
                     tlen
             );
##### add alignments to mmseqs results table?
if( $alnSeqs eq "T" ) {
    push(@mmfields,"qseq","tseq");
}
my $mmfields = join(",",@mmfields);

my $blastOptions
    = qq( -evalue $maxEvalue )
    . qq( -seg yes )
    . qq( -soft_masking true )
    . qq( -comp_based_stats 0 )
    . qq( -parse_deflines )
    . qq( -num_threads $cpus )
    . qq( -outfmt '$pwTbl' );

my $diamondOptions
    = qq( --evalue $maxEvalue )
    . qq( --masking 0 )
    . qq( --comp-based-stats 0 )
    . qq( --threads $cpus )
    . qq( --tmpdir $tempFolder )
    . qq( --quiet )
    . qq( --outfmt $pwTbl );
if( $dmdmode =~ m{sensitive} ) {
    $diamondOptions = qq( --$dmdmode ) . $diamondOptions;
}

my $mmseqsOptions
    = qq( -e $maxEvalue )
    . qq( --comp-bias-corr 0 )
    . qq( --threads $cpus )
    . qq( --format-output "$mmfields" );

#############################################################
################# running the RBH comparisons ###############
#############################################################
my $maxKeep  = 30;
my $kept     = 0;
if( $allvsall > 0 ) {
    ##### by running these in reverse we ensure that
    ##### each database is formatted only once
  TARGETAVA:
    while( my $faaTarget = pop @queries ) {
        my $cntTargets = checkFasta("$faaTarget");
        my $maxAlns
            = $topMatch >= $minTop ? $topMatch
            : sprintf("%.0f",($cntTargets/$matchRatio));
        if( $cntTargets < 1 ) {
            print "  jumping empty file: $faaTarget";
            next TARGETAVA;
        }
        else {
            for my $faaQuery ( @queries,$faaTarget ) {
                runPairWise("$faaQuery","$faaTarget","$maxAlns");
                if( $runRBH eq 'T' ) {
                    buildNrunRBH("$faaQuery","$faaTarget");
                }
            }
            $kept++;
            if( $kept >= $maxKeep  ) {
                system "rm -r $tempFolder/* &>/dev/null";
                $kept = 0;
            }
        }
    }
}
else {
    ##### by running these target-wise to try and ensure that
    ##### each database is formatted only once
  TARGETQVA:
    while( my $faaTarget = shift @against ) {
        my $cntTargets = checkFasta("$faaTarget");
        my $maxAlns
            = $topMatch >= $minTop ? $topMatch
            : sprintf("%.0f",($cntTargets/$matchRatio));
        if( $cntTargets < 1 ) {
            print "  jumping empty file: $faaTarget";
            next TARGETQVA;
        }
        else {
            for my $faaQuery ( @queries ) {
                runPairWise("$faaQuery","$faaTarget","$maxAlns");
                if( $runRBH eq 'T' ) {
                    buildNrunRBH("$faaQuery","$faaTarget");
                }
            }
            $kept++;
            if( $kept >= $maxKeep  ) {
                system "rm -r $tempFolder/* &>/dev/null";
                $kept = 0;
            }
        }
    }
}
##############################################################
######### finish
##############################################################
if( -d "$tempFolder" ) {
    print "   cleaning up\n";
    system("rm -r $tempFolder");
}
print "      Done with $ownName!\n\n";

#########################################################################
############################ subroutines ################################
#########################################################################

sub runPairWise {
    my($faaQuery,$faaTarget,$maxAlns) = @_;
    my $queryGnm   = nakedName($faaQuery);
    my $targetGnm = nakedName($faaTarget);
    my $pwDB       = $tempFolder . "/$targetGnm";
    my $alnFile    = $pwDir . "/$queryGnm/$targetGnm." . $pwProg . ".bz2";
    ### check if we have pairwise comparison results
    ### decide if this should run
    if( -f "$alnFile" && $keepold =~ m{T|R}) {
        print "we already have an alignment file:\n $alnFile\n";
    }
    else {
        print "comparing sequences with $pwProg\n";
        if( $pwProg eq "blastp" ) {
            runBlastp("$faaQuery","$faaTarget","$alnFile","$maxAlns");
        }
        elsif( $pwProg eq "diamond" ) {
            runDiamond("$faaQuery","$faaTarget","$alnFile","$maxAlns");
        }
        elsif( $pwProg eq "mmseqs" ) {
            runMMseqs("$faaQuery","$faaTarget","$alnFile","$maxAlns");
        }
        elsif( $pwProg eq "lastal" ) {
            runLastal("$faaQuery","$faaTarget","$alnFile","$maxAlns");
        }
        else {
            print "no $pwProg program\n\n";
            signalHandler();
        }
    }
}

sub buildNrunRBH {
    my($faaQuery,$faaTarget) = @_;
    my $queryGnm  = nakedName($faaQuery);
    my $targetGnm = nakedName($faaTarget);
    my $pwDB      = $tempFolder . "/$targetGnm";
    my $alnFile   = $pwDir . "/$queryGnm/$targetGnm." . $pwProg . ".bz2";
    #########################################################################
    ######### reciprocal best hits
    #########################################################################
    my $qRBHdir  = "$rbhDir/$queryGnm";
    my $tRBHdir  = "$rbhDir/$targetGnm";
    my $qRBHfile = "$qRBHdir/$targetGnm.rbh.bz2";
    my $tRBHfile = "$tRBHdir/$queryGnm.rbh.bz2";
    if( -f "$qRBHfile" && -f "$tRBHfile" && $keepold =~ m{T} ) {
        print join("\n ","there is RBH files already:"
                       ,$qRBHfile,$tRBHfile),"\n";
    }
    else {
        print "gathering reciprocal best hits\n";
        ## learn top homologs
        print "   learning top homologs\n";
        ##### first learn lines and their scores from query genome:
        my %pair_query  = (); # obvious
        my %pair_target = (); # obvious
        my %pair_bitsc  = (); # capture bit score
        my %pair_Evalue = (); # capture E-value
        my %pair_line   = (); # capture line
        my %ident_bitsc = (); # capture bit score for identical proteins
        my %ident_qline = (); # capture line for identical proteins
        my %ident_tline = (); # capture line for identical proteins
        ##### open original homologs
        open( my $IN,"-|","bzip2 -qdc $alnFile" );
      HOMOLOGLINE:
        while(<$IN>) {
            next HOMOLOGLINE if( m{^#} );
            chomp;
            my(
                $query,$target,
                $evalue,$bitscore,$pident,
                $qstart,$qend,$qlen,
                $sstart,$send,$slen,
                $qalnseq,$salnseq
            ) = split;
            $query   = cleanID("$query");
            $target = cleanID("$target");
            ## percent coverages
            my $qcov = calcCoverage($qstart,$qend,$qlen);
            my $scov = calcCoverage($sstart,$send,$slen);
            if( ( $qcov >= $minCov ) || ( $scov >= $minCov ) ) {
                my $stats = join("\t",
                                 $evalue,$bitscore,
                                 $qstart,$qend,$qcov,
                                 $sstart,$send,$scov
                             );
                if( ( $qlen eq $slen )
                        && ( $pident >= 100 )
                        && ( $qend   eq $qlen )
                        && ( $qstart eq $sstart )
                        && ( $qend   eq $send )
                        && ( $bitscore > $ident_bitsc{"$query"} ) ) {
                    $ident_bitsc{"$query"} = $bitscore;
                    $ident_qline{"$query"}
                        = join("\t",$query,$target,$stats);
                    $ident_tline{"$target"}
                        = join("\t",$target,$query,$stats);
                }
                my $pair = join("::",$query,$target);
                if( $bitscore > $pair_bitsc{"$pair"} ) {
                    $pair_line{"$pair"}   = join("\t",$query,$target,$stats);
                    $pair_query{"$pair"}  = $query;
                    $pair_target{"$pair"} = $target;
                    $pair_bitsc{"$pair"}  = $bitscore;
                    $pair_Evalue{"$pair"} = $evalue;
                }
            }
        }
        close($IN);
        #### now sort, first bit-score (highest to lowest),
        #### then E value (lowest to highest), then query, then target
        my @query_lines = ();
        my @query_pairs
            = sort {
                $pair_bitsc{"$b"} <=> $pair_bitsc{"$a"}
                    || $pair_Evalue{"$a"} <=> $pair_Evalue{"$b"}
                    || $pair_query{"$a"}  <=> $pair_query{"$b"}
                    || $pair_query{"$a"}  cmp $pair_query{"$b"}
                    || $pair_target{"$a"} <=> $pair_target{"$b"}
                    || $pair_target{"$a"} cmp $pair_target{"$b"}
                } keys %pair_query;
        for my $pair ( @query_pairs ) {
            push(@query_lines,$pair_line{"$pair"});
        }
        #### now sort reciprocals
        my @reciprocal_lines = ();
        my @pretend_reciprocal
            = sort {
                $pair_bitsc{"$b"} <=> $pair_bitsc{"$a"}
                    || $pair_Evalue{"$a"} <=> $pair_Evalue{"$b"}
                    || $pair_target{"$a"} <=> $pair_target{"$b"}
                    || $pair_target{"$a"} cmp $pair_target{"$b"}
                    || $pair_query{"$a"}  <=> $pair_query{"$b"}
                    || $pair_query{"$a"}  cmp $pair_query{"$b"}
                } keys %pair_query;
        for my $pair ( @pretend_reciprocal ) {
            my (
                $query,$target,
                $evalue,$bit_score,
                $q_start,$q_end,$qcov,
                $s_start,$s_end,$scov
            ) = split(/\t/,$pair_line{"$pair"});
            my $oppLine = join("\t",
                               $target,$query,
                               $evalue,$bit_score,
                               $s_start,$s_end,$scov,
                               $q_start,$q_end,$qcov
                           );
            push(@reciprocal_lines,$oppLine);
        }
        ###### now orthology
        ### directories for RBHs:
        my $tmpQuery  = $tempFolder . "/$queryGnm.$targetGnm.rbh.bz2";
        my $tmpTarget = $tempFolder . "/$targetGnm.$queryGnm.rbh.bz2";
        print "   extracting orthologs\n";
        ###### making this into a subroutine to allow working both ways:
        produceRBH($tmpQuery,\@query_lines,\@reciprocal_lines,
                   \%ident_qline,\%ident_tline);
        if( -s "$tmpQuery" ) {
            for my $checkDir ( "$rbhDir","$qRBHdir" ) {
                unless( -d "$checkDir" ) {
                    mkdir("$checkDir");
                }
            }
            system( "mv $tmpQuery $qRBHfile 2>/dev/null" );
        }
        produceRBH($tmpTarget,\@reciprocal_lines,\@query_lines,
                   \%ident_tline,\%ident_qline);
        if( -s "$tmpTarget" ) {
            for my $checkDir ( "$rbhDir","$tRBHdir" ) {
                unless( -d "$checkDir" ) {
                    mkdir("$checkDir");
                }
            }
            system( "mv $tmpTarget $tRBHfile 2>/dev/null" );
        }
    }
}

sub produceRBH {
    my ($outFile,$rqLines,$rsLines,$ridentLnQ,$ridentLnS) = @_;
    my $tmpFile = $outFile . ".tmp";
    my %query_best_hits = ();
    my %query_bitsc     = ();
    my %query_E_value   = ();
    my %init_target     = ();
    my %end_target      = ();
    my %print_line      = ();
  QLINE:
    for my $line ( @{$rqLines} ) {
        my ( $query,$target,$evalue,$bit_score,
             $q_start,$q_end,$qcov,
             $s_start,$s_end,$scov
         ) = split(/\s+/,$line);
        next QLINE if( exists $ridentLnQ->{"$query"} ); ### jump identicals
        my $print_line = join("\t",
                              $query,$target,
                              $evalue,$bit_score,
                              $q_start,$q_end,$qcov,
                              $s_start,$s_end,$scov
                          );
        if( length($query_best_hits{"$query"}) > 0 ) {
            if( $query_bitsc{"$query"} == $bit_score ) {
                $query_best_hits{"$query"}    .= "," . $target;
                $print_line{"$query.$target"}  = $print_line;
                $init_target{"$query.$target"} = $s_start;
                $end_target{"$query.$target"}  = $s_end;
            }
        }
        else {
            $query_best_hits{"$query"}     = $target;
            $print_line{"$query.$target"}  = $print_line;
            $init_target{"$query.$target"} = $s_start;
            $end_target{"$query.$target"}  = $s_end;
            $query_bitsc{"$query"}         = $bit_score;
            $query_E_value{"$query"}       = $evalue;
        }
    }
    #####
    my %target_best_hits = ();
    my %target_bitsc     = ();
    my %target_E_value   = ();
  SLINE:
    for my $line ( @{$rsLines} ) {
        my ( $query,$target,$evalue,$bit_score,
             $q_start,$q_end,$qcov,
             $s_start,$s_end,$scov
         ) = split(/\s+/,$line);
        next SLINE if( exists $ridentLnS->{"$query"} ); ### jump identicals
        if( length $target_best_hits{"$query"} > 0 ) {
            if( $target_bitsc{"$query"} == $bit_score ) {
                $target_best_hits{"$query"} .= "," . $target;
            }
        }
        else {
            $target_best_hits{"$query"} = $target;
            $target_bitsc{"$query"}     = $bit_score;
            $target_E_value{"$query"}   = $evalue;
        }
    }
    #####
    my %orth           = ();
    my %share_this_hit = ();
    ##### sort the queries by bit score to automatically sort the
    ##### best hits by their bit scores (I guess)
    my @queries = sort {
        $query_bitsc{"$b"} <=> $query_bitsc{"$a"}
            || $query_E_value{"$a"} <=> $query_E_value{"$b"}
            || $a <=> $b
            || $a cmp $b
        } keys(%query_best_hits);
    ##### now do the finding of orthologs
    for my $query ( @queries ) {
        for my $query_best_hit ( split(/\,/,$query_best_hits{"$query"}) ) {
            for my $target_best_hit (
                split(/\,/,$target_best_hits{"$query_best_hit"})
            ) {
                ## Reciprocal best hits:
                if( $target_best_hit eq $query ) {
                    $orth{"$query"} .= ",$query_best_hit";
                    $orth{"$query"} =~ s/^\,+//;
                }
                else {
                    $share_this_hit{"$query_best_hit"} .= ",$query";
                    $share_this_hit{"$query_best_hit"} =~ s/^\,+//;
                }
            }
        }
    }
    ##### count printed lines (to ensure content)
    my $printedOrths = 0;
    open( my $ORTHS,"|-","bzip2 -9 >$tmpFile" );
    print {$ORTHS} "#",join("\t",@rbhHeading),"\n";
    ###### first print identicals:
    for my $query ( sort { $a <=> $b || $a cmp $b } keys %{$ridentLnQ} ) {
        print {$ORTHS} join("\t",$ridentLnQ->{"$query"},"Identical"),"\n";
        $printedOrths++;
    }
    ######
    ###### now reciprocal best hits, and fusions:
    for my $query ( sort { $a <=> $b || $a cmp $b } keys %orth ) {
      ORTH:
        for my $orth ( split(/,+/,$orth{"$query"}) ) {
            #### print reciprocal best hit ortholog
            print {$ORTHS} $print_line{"$query.$orth"},"\tRBH\n";
            $printedOrths++;
            ### find fusions from query to target by testing
            ### other queries that find the same target as their top
            ### hit, but not reciprocal
            my @share_this_hit = split(/\,/,$share_this_hit{"$orth"});
            my $shared = @share_this_hit;
            next ORTH unless( $shared > 0 );
            my ($i_query,$f_query)
                = ($init_target{"$query.$orth"},$end_target{"$query.$orth"});
            my $rbh_coords = join(":",$i_query,$f_query);
            my @coords = ("$rbh_coords");
          PARTNER:
            for my $partner ( @share_this_hit ) {
                #print $print_line{"$partner.$orth"},"\n";
                my $q_start = $init_target{"$partner.$orth"};
                my $q_end   = $end_target{"$partner.$orth"};
                ##### allowing a bit of overlap for fusions:
                my $q_max_o
                    = sprintf("%.0f",$maxOverlap * ($q_end - $q_start));
              COORDS:
                for my $coords ( @coords ) {
                    my($cover_start,$cover_end) = split(/:/,$coords);
                    ##### allowing a bit of overlap:
                    my $cover_max_o
                        = sprintf("%.0f",
                                  $maxOverlap*($cover_end - $cover_start));
                    my $min_overlap
                        = $cover_max_o > $q_max_o ? (-1 * $q_max_o)
                        : (-1 * $cover_max_o);
                    my $overlap1 = $q_start - ( $cover_end + 1);
                    my $overlap2 = $cover_start - ( $q_end + 1);
                    #print "OVERLAPS: $overlap1, $overlap2, $min_overlap\n";
                    if( $overlap1 < $min_overlap
                            && $overlap2 < $min_overlap ) {
                        next PARTNER;
                    }
                    #print "PASSED\n";
                }
                my $new_coords = join(":",$q_start,$q_end);
                push(@coords,$new_coords);
                if( $q_start < $i_query ) {
                    print {$ORTHS}
                        $print_line{"$partner.$orth"}
                        ,"\tLEFT of " . $query . "\n";
                    $printedOrths++;
                }
                else {
                    print {$ORTHS}
                        $print_line{"$partner.$orth"}
                        ,"\tRIGHT of " . $query . "\n";
                    $printedOrths++;
                }
            } # for partner
        } # for orth
    } # for query
    close($ORTHS);
    ### rename resulting file, then check if it has content
    ### erase otherwise
    rename($tmpFile,$outFile);
    if( -z "$outFile" || $printedOrths == 0 ) {
        print "   found no orthologs\n";
        unlink("$outFile");
    }
}

sub nakedName {
    my $inName  = $_[0];
    my $outName = $inName;
    $outName =~ s{^\S+/}{};
    $outName =~ s{\.(Z|gz|bz2)$}{};
    $outName =~ s{\.faa\S*}{};
    $outName =~ s{\.fasta\S*}{};
    if( length($outName) > 0 ) {
        return("$outName");
    }
    else {
        print "there's a problem with $inName -> $outName\n";
        return();
    }
}

sub nameDB {
    my $inName  = $_[0];
    my $outName = nakedName("$inName");
    if( length($outName) > 0 ) {
        return("$tempFolder/$outName");
    }
    else {
        print "there's a problem with $inName -> $outName\n";
        return();
    }
}

sub signalHandler {
    if( length($tempFolder) > 1 && -d "$tempFolder" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $tempFolder";
    }
    else {
        print "\n\tquitting $ownName\n";
    }
    die  "\tdone!\n\n";
}

sub how2open {
    my $rootName = $_[0];
    $rootName =~ s{\.(gz|bz2|Z)$}{};
    my $fullName
        = ( -f "$rootName.gz" )  ? "$rootName.gz"
        : ( -f "$rootName.Z" )   ? "$rootName.Z"
        : ( -f "$rootName.bz2" ) ? "$rootName.bz2"
        : ( -f "$rootName" )     ? "$rootName"
        : "none";
    my $catProg
        = $fullName =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : $fullName =~ m{\.bz2$}    ? "bzip2 -qdc"
        : "cat";
    if( $fullName eq "none" ) {
        return();
    }
    else {
        return("$catProg","$fullName");
    }
}

sub inflateFile {
    my $origfile = $_[0];
    my( $opener ) = how2open("$origfile");
    my $rootName = nakedName("$origfile");
    my $inflated = "$tempFolder/$rootName.faa";
    unless( -f "$inflated" ) {
        system("$opener $origfile > $inflated");
    }
    return("$inflated");
}

sub formatDB {
    my($file,$dbfile,$dbType) = @_;
    my($opener) = how2open("$file");
    if( $dbType eq "blastp" ) {
        if( -f "$dbfile.psq" ) {
            print "  the blast DB file is already there\n";
        }
        else {
            print "producing blastDB: $dbfile\n";
            my $mkblastdb
                = qq($opener $file |)
                . qq( makeblastdb )
                . qq( -dbtype prot )
                . qq( -title $dbfile )
                . qq( -parse_seqids )
                . qq( -out $dbfile );
            system("$mkblastdb 1>/dev/null");
        }
    }
    elsif( $dbType eq "diamond" ) {
        if( -f "$dbfile.dmnd" ) {
            print "  the diamond DB file is already there\n";
        }
        else {
            print "producing diamondDB: $dbfile\n";
            my $mkdiamondDB
                = qq($opener $file |)
                . qq( diamond makedb -d $dbfile --quiet --threads $cpus );
            system("$mkdiamondDB &>/dev/null");
        }
    }
    elsif( $dbType eq "mmseqs" ) {
        if( -f "$dbfile.index" ) {
            print "  the mmseqs DB file is already there\n";
        }
        else {
            print "producing mmseqsDB: $dbfile\n";
            my $mkmmseqsDB
                = qq( mmseqs createdb $file $dbfile );
            system("$mkmmseqsDB &>/dev/null");
            #seems like the idx files make mmseqs less sensitive
            #print "producing $dbfile mmseqs index\n";
            #my $mkIndex = qq( mmseqs createindex $dbfile $tempFolder );
            #system("$mkIndex &>/dev/null");
        }
    }
    elsif( $dbType eq "lastal" ) {
        if( -f "$dbfile.suf" ) {
            print "  the lastal DB file is already there\n";
        }
        else {
            print "producing lastalDB: $dbfile\n";
            my $mklastalDB
                = qq($opener $file | segmasker -outfmt fasta |)
                . qq( lastdb -p -c $dbfile );
            system("$mklastalDB &>/dev/null");
        }
    }
    else {
        print "I need a database format [$matchProg]\n\n";
    }
}

sub calcCoverage {
    my($start,$end,$ln) = @_;
    if( $ln < 1 ) {
        return();
    }
    else {
        my $coverage = 100 * ( ( $end - $start + 1 ) / $ln );
        my $rounded  = sprintf( "%.1f", $coverage );
        return($rounded);
    }
}

sub runBlastp {
    my($queryFile,$targetFile,$alnFile,$maxAlns) = @_;
    my $dbfile = nameDB("$targetFile");
    formatDB("$targetFile","$dbfile","$pwProg");
    my $nkQuery   = nakedName("$queryFile");
    my $nkTarget = nakedName("$targetFile");
    my $tmpOut = "$tempFolder/$nkQuery.$nkTarget.$pwProg.bz2";
    print "running blastp:\n   ",$nkQuery," vs ",$nkTarget,"\n";
    my( $opener,$file ) = how2open($queryFile);
    my $blcommand
        = qq($opener $file | blastp -db $dbfile )
        . qq(-max_target_seqs $maxAlns $blastOptions);
    my $lineCount = 0;
    open( my $PWBLAST,"|-","bzip2 -9 > $tmpOut" );
    print {$PWBLAST} "# ",join("\t",@pwHeading),"\n";
    for my $blastLine ( qx($blcommand) ) {
        $lineCount++;
        print {$PWBLAST} $blastLine;
    }
    close($PWBLAST);
    if( $lineCount > 0 ) {
        unless( -d "$pwDir" ) {
            mkdir("$pwDir");
        }
        unless( -d "$pwDir/$nkQuery" ) {
            mkdir("$pwDir/$nkQuery");
        }
        system qq(mv $tmpOut $alnFile &>/dev/null);
        print "   ",$nkQuery," vs ",$nkTarget," done\n";
    }
    else {
        print "   blastp failed (empty file)\n";
        signalHandler();
    }
}

sub runDiamond {
    my($queryFile,$targetFile,$alnFile,$maxAlns) = @_;
    my $dbfile = nameDB("$targetFile");
    formatDB("$targetFile","$dbfile","$pwProg");
    my $nkQuery   = nakedName("$queryFile");
    my $nkTarget = nakedName("$targetFile");
    my $tmpOut = "$tempFolder/$nkQuery.$nkTarget.$pwProg.bz2";
    print "running diamond:\n   ",$nkQuery," vs ",$nkTarget,"\n";
    my( $opener,$file ) = how2open("$queryFile");
    my $dmcommand
        = qq($opener $file | diamond blastp --db $dbfile )
        . qq(--max-target-seqs $maxAlns $diamondOptions);
    my $lineCount = 0;
    open( my $PWDMD,"|-","bzip2 -9 > $tmpOut" );
    print {$PWDMD} "# ",join("\t",@pwHeading),"\n";
    for my $dmLine ( qx($dmcommand 2>/dev/null) ) {
        $lineCount++;
        print {$PWDMD} $dmLine;
    }
    close($PWDMD);
    if( $lineCount > 0 ) {
        unless( -d "$pwDir" ) {
            mkdir("$pwDir");
        }
        unless( -d "$pwDir/$nkQuery" ) {
            mkdir("$pwDir/$nkQuery");
        }
        system qq(mv $tmpOut $alnFile &>/dev/null);
        print "   ",$nkQuery," vs ",$nkTarget," done\n";
    }
    else {
        print "   diamond failed (empty file)\n";
        signalHandler();
    }
}

sub runMMseqs {
    my($queryFile,$targetFile,$alnFile,$maxAlns) = @_;
    my $dbfile = nameDB("$targetFile");
    formatDB("$targetFile","$dbfile","$pwProg");
    my $nkQuery   = nakedName("$queryFile");
    my $nkTarget = nakedName("$targetFile");
    print "running mmseqs:\n   ",$nkQuery," vs ",$nkTarget,"\n";
    my $tmpOut = "$tempFolder/$nkQuery.$nkTarget.tmp";
    my $mmcommand
        = qq( mmseqs easy-search $queryFile $dbfile $tmpOut $tempFolder )
        . qq( --max-accept $maxAlns $mmseqsOptions);
    my $mmLog = qx($mmcommand 2>&1);
    my $lineCount = 0;
    open( my $MMOUT,"<","$tmpOut" );
    open( my $PWMM,"|-","bzip2 -9 > $tmpOut.bz2" );
    print {$PWMM} "# ",join("\t",@pwHeading),"\n";
    while(<$MMOUT>) {
        $lineCount++;
        my( $query,$target,
            $evalue,$bits,$pident,
            $qstart,$qend,$qlen,
            $tstart,$tend,$tlen,
            $qseq,$tseq
        ) = split;
        ##### fix percent identity
        $pident *= 100;
        my @putBack = (
            $query,$target,
            $evalue,$bits,$pident,
            $qstart,$qend,$qlen,
            $tstart,$tend,$tlen
        );
        if( length($qseq) > 0 && length($tseq) > 0 ) {
            push(@putBack,$qseq,$tseq);
        }
        print {$PWMM} join("\t",@putBack),"\n";
    }
    close($MMOUT);
    close($PWMM);
    if( $lineCount > 0 ) {
        unless( -d "$pwDir" ) {
            mkdir("$pwDir");
        }
        unless( -d "$pwDir/$nkQuery" ) {
            mkdir("$pwDir/$nkQuery");
        }
        system qq(mv $tmpOut.bz2 $alnFile &>/dev/null);
        print "   ",$nkQuery," vs ",$nkTarget," done\n";
    }
    else {
        print "   mmseqs failed (empty file)\n";
        signalHandler();
    }
}

sub runLastal {
    my($queryFile,$targetFile,$alnFile,$maxAlns) = @_;
    my $dbfile = nameDB("$targetFile");
    formatDB("$targetFile","$dbfile","$pwProg");
    my $nkQuery   = nakedName("$queryFile");
    my $nkTarget = nakedName("$targetFile");
    my $tmpOut = "$tempFolder/$nkQuery.$nkTarget.$pwProg.bz2";
    print "running lastal:\n   ",$nkQuery," vs ",$nkTarget,"\n";
    my( $opener,$file ) = how2open("$queryFile");
    my $lastcommand
        = qq($opener $file | lastal -f BlastTab+ -P $cpus -N $maxAlns $dbfile);
    my $lineCount = 0;
    open( my $PWLAST,"|-","bzip2 -9 > $tmpOut" );
    print {$PWLAST} "# ",join("\t",@pwHeading),"\n";
  LASTLINE:
    for my $lastLine ( qx($lastcommand 2>/dev/null) ) {
        next LASTLINE if( $lastLine =~ m{^#} );
        next LASTLINE if( $lastLine =~ m{^\s*\n} );
        chomp;
        my( $query,$target,$pident,
            $alignmentLength,$mismatches,$gapopens,
            $qstart,$qend,
            $tstart,$tend,
            $evalue,$bits,
            $qlen,$tlen,$raw
        ) = split(/\t/,$lastLine);
        #next LASTLINE if ( $evalue > $maxEvalue );
        $lineCount++;
        print {$PWLAST}
            join("\t",
                 $query,$target,
                 $evalue,$bits,$pident,
                 $qstart,$qend,$qlen,
                 $tstart,$tend,$tlen),"\n";
    }
    close($PWLAST);
    if( $lineCount > 0 ) {
        unless( -d "$pwDir" ) {
            mkdir("$pwDir");
        }
        unless( -d "$pwDir/$nkQuery" ) {
            mkdir("$pwDir/$nkQuery");
        }
        system qq(mv $tmpOut $alnFile &>/dev/null);
        print "   ",$nkQuery," vs ",$nkTarget," done\n";
    }
    else {
        print "   lastal failed (empty file)\n";
        signalHandler();
    }
}

sub cleanID {
    my $toClean = $_[0];
    $toClean =~ s{^\S+?\|(\S+)}{$1};
    $toClean =~ s{\|$}{};
    return("$toClean");
}

sub findFaaFiles {
    my($testDir,$label) = @_;
    print "finding $label faa|fasta files in $testDir\n";
    my @foundFiles = ();
    #### now repopulate arrays:
    opendir( my $FAADIR,"$testDir" );
    my @faaFiles = grep { m{\.(faa|fasta)} } readdir $FAADIR;
    for my $faaFile ( @faaFiles ) {
        push( @foundFiles,"$testDir/$faaFile" );
    }
    my $cntFound = @foundFiles;
    if( $cntFound > 0 ) {
        print " will work with $cntFound $label files\n";
        return(@foundFiles);
    }
    else {
        die " no faa/fasta files in $label directory ($testDir)\n\n";
    }
}

sub checkFasta {
    my $checkFaa = $_[0];
    my( $opener,$file ) = how2open("$checkFaa");
    my $cntSeqs = 0;
    open( my $CFAA,"-|","$opener $file" );
    while(<$CFAA>) {
        if( m{^>} ) {
            $cntSeqs++;
        }
    }
    close($CFAA);
    return($cntSeqs);
}

sub findWorkFiles {
    my ($item,@testFiles) = @_;
    my $toFind = @testFiles;
    print "checking for $toFind $item files\n";
    my @found = ();
    for my $file ( sort @testFiles ) {
        if( -f "$file" ) {
            push(@found,$file);
        }
        else {
            print "did not find $file\n";
        }
    }
    my $foundFiles = @found;
    my $missing = $toFind - $foundFiles;
    if( $missing > 0 ) {
        if( $missing == 1 ) {
            die " I could not find $missing $item file\n\n";
        }
        else {
            die " I could not find $missing $item files\n\n";
        }
    }
    else {
        return(@found);
    }
}

sub checkSoftware {
    my %available = ();
    my %missing   = ();
    for my $prog ( @pwProgs ) {
        my $path2p = qx(which $prog 2>&1);
        chomp($path2p);
        if( $path2p =~ m{^\S+/$prog$} ) {
            $available{"$prog"}++;
        }
        else {
            $missing{"$prog"}++;
        }
    }
    my @available  = sort keys %available;
    my $availCount = @available;
    my @missing    = sort keys %missing;
    my $missCount  = @missing;
    my $missMsg    = '';
    if( $missCount > 0 ) {
        my $num = 1;
        for my $missing ( @missing ) {
            $missMsg
                .= $num . " " . $missing
                . " available at:\n";
            $missMsg .= $url{"$missing"} . "\n\n";
            $num++;
        }
    }
    if( $availCount < 1 ) {
        print qq(\nThere's no software to find RBHs:\n\n),$missMsg;
        exit;
    }
    else {
        return(\@available,\%missing,"$missMsg");
    }
}

sub podhelp {
    my $extraMessage = $_[0];
    open( my $PIPE,"|-","less -wiseMR" );
    if( length "$extraMessage" > 2 ) {
        print {$PIPE} $extraMessage,"\n\n";
    }
    $parser->output_fh($PIPE);
    $parser->parse_string_document($podUsage);
    exit;
}
