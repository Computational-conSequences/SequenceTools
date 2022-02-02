#!/usr/bin/perl
use strict;
use List::Util qw(sum);
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my $ownName = $0;
$ownName =~ s{.*/}{};

my %url = (
    "ANIx" => 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/',
    "ANIp" => 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/',
    "ANIb" => 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/',
    "ANIm" => 'https://github.com/mummer4/mummer/releases/',
    "ANIf" => 'https://github.com/ParBLiSS/FastANI/releases/',
    "ANIl" => 'http://last.cbrc.jp/',
    "ANIu" => 'https://www.drive5.com/usearch/',
);
my @methods = sort keys %url;
my $methods = join("|",@methods);

my %pwProg = (
    "ANIx" => 'blastp',
    "ANIp" => 'blastp',
    "ANIb" => 'blastall',
    "ANIm" => 'nucmer',
    "ANIf" => 'fastANI',
    "ANIl" => 'lastal',
    "ANIu" => 'usearch',
);

######## check for software to calculate ANI:
my %available = ();
my %missing   = ();
for my $testM ( sort keys %pwProg ) {
    my $prog = $pwProg{"$testM"};
    my $path2p = qx(which $prog 2>&1);
    chomp($path2p);
    if( $path2p =~ m{^\S+/$prog$} ) {
        $available{"$testM"}++;
    }
    else {
        $missing{"$testM"}++;
    }
}
my @available  = sort keys %available;
my $available  = join("|",@available);
my $availCount = @available;
my @missing    = sort keys %missing;
my $missing    = join("|",@missing);
my $missCount  = @missing;
my $missMsg    = '';
if( $missCount > 0 ) {
    my $num = 1;
    for my $missing ( @missing ) {
        $missMsg
            .= $num . ". " . $missing . ": " . $pwProg{"$missing"}
            . " available at:\n";
        $missMsg .= $url{"$missing"} . "\n\n";
        $num++;
    }
}
if( $availCount < 1 ) {
    print qq(\nThere's no software to calculate ANI:\n\n),$missMsg;
    exit;
}

my %needCuts = (
    "ANIb" => 1,
    "ANIu" => 1,
    "ANIl" => 1,
    "ANIp" => 1,
);

my $cpuNumber
    = qx(getconf _NPROCESSORS_ONLN 2>/dev/null)
    =~ m{\s*(\d+)\s*}                  ? $1
    : qx(sysctl -a 2>/dev/null | grep 'cpu.thread_count')
    =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
    =~ m{\.max-threads\s+=\s+(\d+)}    ? $1
    : 2;

my $defCPUs
    = $cpuNumber >= 2 ? 2
    : 1;

my $minLn      = 1000;
my $maxLn      = 3000;
my $defaultLn  = 1020;
my $defaultANI
    = exists $available{"ANIf"} ? "ANIf"
    : exists $available{"ANIx"} ? "ANIx"
    : exists $available{"ANIm"} ? "ANIm"
    : exists $available{"ANIp"} ? "ANIp"
    : $available[0];
####### setting program arguments:
my @queries    = ();
my @against    = ();
my $fnaDir     = '';
my $method     = $defaultANI;
my $cutLn      = $defaultLn;
my $resultsDir = 'ResultsANI';
my $cpus       = $defCPUs;
my $keepold    = 'T';
####### when working on all vs all we use a
####### number of files for batch work:
my $batch      = 200;

my $helpMsg
    = qq(about:\n)
    . qq(  This program calculates Average Nucleotide Identity\n)
    . qq(  (ANI) by several methods\n\n)
    . qq(usage: $ownName -q <fnafile> -t <fnafile> [options]\n\n)
    . qq(options:\n)
    . qq(   -q query file in fasta format, or directory with query\n)
    . qq(       files in fasta format, required\n)
    . qq(       (files can be compressed with gzip or bzip2)\n)
    . qq(   -t target file in fasta format, or directory with target\n)
    . qq(       files in fasta format, required\n)
    . qq(       (files can be compressed with gzip or bzip2)\n)
    . qq(   -d directory with fna files for all-vs-all, no default,\n)
    . qq(       if set, -i and -t will be ignored\n)
    . qq(   -m method [$methods], default $defaultANI\n)
    . qq(       ANIx: our blast+ reinterpretation (megablast)\n)
    . qq(       ANIp: our blast+ implementation of original method (slow)\n)
    . qq(       ANIb: original method using legacy blast (slowest)\n)
    . qq(             https://www.doi.org/10.1099/ijs.0.64483-0\n)
    . qq(       ANIm: mummer method (using nucmer)\n)
    . qq(             https://www.doi.org/10.1073/pnas.0906412106\n)
    . qq(       ANIf: fastANI method\n)
    . qq(             https://www.doi.org/10.1038/s41467-018-07641-9\n)
    . qq(       ANIl: our lastal implementation\n)
    . qq(       ANIu: our ublast implementation\n)
    . qq(   -c fragment length (not used in ANIx, or ANIm),\n)
    . qq(       minimum $minLn, maximum $maxLn default $defaultLn\n)
    . qq(   -o output folder, default ResultsANI\n)
    . qq(   -k keep prior result [T|F], default 'T'. If 'T' prior results\n)
    . qq(       will be taken from appropriate files in the output folder\n)
    . qq(   -x number of cpus to use, maximum $cpuNumber, default $defCPUs\n)
    #. qq(       using '-x X' with ANIf will distribute genomes into lists of\n)
    #. qq(       $batch genomes and run fastANI with the maximum number of cpus.\n)
    . qq(\n);
if( $missCount > 0 ) {
    $helpMsg
        .= qq(Warning:\n)
        .  qq(  The software for the following methods is missing:\n\n)
        .  $missMsg;
}

my $options = GetOptions(
    "q=s{,}" => \@queries,
    "t=s{,}" => \@against,
    "d=s"    => \$fnaDir,
    "m=s"    => \$method,
    "c=i"    => \$cutLn,
    "o=s"    => \$resultsDir,
    "k=s"    => \$keepold,
    "x=i"    => \$cpus,
) or die "$helpMsg";

### check if working with all-vs-all directory
my $allvsall = length($fnaDir) > 0 ? 1 : 0;
###### check files acccordingly
if( $allvsall == 0 ) {
    if( !$queries[0] || !$against[0] ) {
        print $helpMsg;
        exit;
    }
    else {
        ####### check queries
        my %qseen = ();
        @queries = sort grep { not $qseen{"$_"}++ } @queries;
        my $qtoFind = @queries;
        if( -d $queries[0] && $qtoFind == 1 ) {
            @queries = findFnaFiles("$queries[0]","query");
        }
        else {
            @queries = findWorkFiles("query",@queries);
        }
        ###### check targets
        my %tseen = ();
        @against = sort grep { not $tseen{"$_"}++ } @against;
        my $stoFind = @against;
        if( -d $against[0] && $stoFind == 1 ) {
            @against = findFnaFiles("$against[0]","target");
        }
        else {
            @against = findWorkFiles("target",@against);
        }
    }
}
else {
    if( -d "$fnaDir" ) {
        if( @queries = findFnaFiles("$fnaDir","all vs all") ) {
            @against = @queries;
        }
    }
    else {
        die "$fnaDir doesn't exist or is not a directory\n\n";
    }
}

my $method = $method =~ m{($methods)} ? $1 : $defaultANI;
if( exists $missing{"$method"} ) {
    my $msgOut
        = qq(The $method method needs missing software:\n)
        . $pwProg{"$method"} . "\n"
        . qq(which is available here:\n)
        . $url{"$method"} . "\n"
        . qq(Otherwise choose an available method [$available]\n);
    die $msgOut;
}
print "using the $method method\n";

my $cutLn = $cutLn >= $minLn && $cutLn <= $maxLn ? $cutLn : $defaultLn;
print "using ${cutLn}bp fragment length\n";
my $keepold = $keepold =~ m{^(T|F)$}i ? uc($1) : 'T';
print "keeping old results if there's any: $keepold\n";
my $cpus    = $cpus > 0 && $cpus <= $cpuNumber ? $cpus : $defCPUs;
print "using $cpus threads\n";
my $maxEvalue = 1e-3;
my $round     = 2;

### directory to save temporary files
my $tempFolder = tempdir("/tmp/$method.XXXXXXXXXXXX");

#################################################################
####### checking for prior results ##############################
#################################################################
print "will save ANI results in the $resultsDir folder\n";
my $qnaked = nakedName($queries[0]);
my $aniout
    = $allvsall > 0 ? "AllvsAll.$method.tbl"
    : "$qnaked.$method.tbl";
$aniout = join("/",$resultsDir,$aniout);
my %prevANI = ();
if( -d "$resultsDir" ) {
    print "ensuring prior result files survival\n";
    opendir( my $PRIORD,"$resultsDir" );
    my @priorFiles = grep { m{\.$method\.} } readdir($PRIORD);
    if( $keepold eq 'T' ) {
        my $saveFile = "prior.$method.backup";
        open( my $PRESERVE,">","$tempFolder/$saveFile" );
        print {$PRESERVE} join("\t","Genome1","Genome2",
                             "G1-G2","G2-G1","ANI"),"\n";
        for my $priorFile (@priorFiles) {
            print "reading previous results ($resultsDir/$priorFile)\n";
            open( my $PREV,"<","$resultsDir/$priorFile" );
          PREVANI:
            while(<$PREV>) {
                next PREVANI if m{^Genome1};
                my ($gnm1,$gnm2,@ani) = split;
                my $pair = join(",",sort($gnm1,$gnm2));
                if( $ani[-1] ne 'NA' && $ani[-1] >= 0.1 ) {
                    unless( exists $prevANI{"$pair"} ) {
                        $prevANI{"$pair"} = join("\t",@ani);
                        print {$PRESERVE}
                            join("\t",sort($gnm1,$gnm2),@ani),"\n";
                    }
                }
            }
            close($PREV);
            unlink("$resultsDir/$priorFile");
        }
        close($PRESERVE);
        system("mv $tempFolder/$saveFile $resultsDir/$saveFile 2>/dev/null");
        my $cntprior = keys %prevANI;
        print "keeping $cntprior prior results\n";
    }
}
else {
    mkdir("$resultsDir");
}

#################################################################
##### Running ANI calculations ##################################
##### trying to make this work with batches of files without
##### getting the main hard drive competely filled
#################################################################
my %batchedFile = ();
my %redundant   = ();
open( my $OUTANI,">","$aniout.tmp" );
print {$OUTANI} join("\t","Genome1","Genome2",
                     "G1-G2","G2-G1","ANI"),"\n";
if( $allvsall > 0 ) {
    while( my $query = shift @queries ) {
        for my $against ( $query,@queries ) {
            my $nakedQ = nakedName("$query");
            my $nakedA = nakedName("$against");
            my $pair = join(",",sort($nakedQ,$nakedA));
            if( exists $prevANI{"$pair"} ) {
                print {$OUTANI}
                    join("\t",sort($nakedQ,$nakedA),$prevANI{"$pair"}),"\n";
            }
            else {
                $batchedFile{"$query"}++;
                $batchedFile{"$against"}++;
                my $batched = keys %batchedFile;
                if( $batched >= $batch ) {
                    %batchedFile = ();
                }
                proceedwQT("$query","$against",$batched);
            }
        }
    }
}
else {
    for my $query ( @queries ) {
      QVTTARGET:
        for my $against ( @against ) {
            my $nakedQ = nakedName("$query");
            my $nakedA = nakedName("$against");
            my $pair = join(",",sort($nakedQ,$nakedA));
            if( exists $redundant{"$pair"} ) {
                next QVTTARGET;
            }
            elsif( exists $prevANI{"$pair"} ) {
                print {$OUTANI}
                    join("\t",sort($nakedQ,$nakedA),$prevANI{"$pair"}),"\n";
            }
            else {
                $batchedFile{"$query"}++;
                $batchedFile{"$against"}++;
                my $batched = keys %batchedFile;
                if( $batched >= $batch ) {
                    %batchedFile = ();
                }
                proceedwQT("$query","$against",$batched);
                $redundant{"$pair"} = 1;
            }
        }
    }
}
close($OUTANI);
rename("$aniout.tmp","$aniout");
if( length($tempFolder) > 1 && -d "$tempFolder" ) {
    print "\tcleaning up ...\n";
    system "rm -r $tempFolder";
}
print "\tdone with $ownName\n\n";

######### first subroutine
sub proceedwQT {
    my( $query,$against,$batched ) = @_;
    print "computing ANI for $query and $against,\n"
        . "  to be saved at: $aniout\n";
    ###### pieces
    if( exists $needCuts{"$method"} ) {
        for my $file ( $query,$against ) {
            my $piecesFile = namePieces("$file");
            print "producing pieces ${cutLn}bp long from\n$file\n";
            cutIntoPieces($file,$piecesFile);
        }
    }
    ###### now run the stuff
    my $ani1
        = $method eq "ANIb" ? calcANIb($query,$against)
        : $method eq "ANIf" ? calcANIf($query,$against)
        : $method eq "ANIu" ? calcANIu($query,$against)
        : $method eq "ANIl" ? calcANIl($query,$against)
        : $method eq "ANIm" ? calcANIm($query,$against)
        : $method eq "ANIp" ? calcANIp($query,$against)
        : calcANIx($query,$against);
    my $ani2
        = $query eq $against ? $ani1
        : $method eq "ANIb"  ? calcANIb($against,$query)
        : $method eq "ANIf"  ? calcANIf($against,$query)
        : $method eq "ANIu"  ? calcANIu($against,$query)
        : $method eq "ANIl"  ? calcANIl($against,$query)
        : $method eq "ANIp"  ? calcANIp($against,$query)
        : $ani1; ### both ANIm and ANIx are reciprocal
        #: $method eq "ANIm" ? calcANIm($against,$query)
        #: calcANIx($against,$query);
    ######## fix results if in trouble:
    if( $ani1 eq '' ) {
        $ani1 = 'NA';
    }
    if( $ani2 eq '' ) {
        $ani2 = 'NA';
    }
    my $ani
        = $ani1 eq 'NA' ? "NA"
        : $ani2 eq 'NA' ? "NA"
        : sprintf("%.${round}f",( ( $ani1 + $ani2 ) / 2 ) );
    my $outLine
        = join("\t",nakedName("$query"),nakedName("$against"),
               $ani1,$ani2,$ani);
    print "ANIs = ",$outLine,"\n";
    print {$OUTANI} $outLine,"\n";
    #### clean up to avoid overloading the hard drive
    if( $batched >= $batch ) {
        system("rm $tempFolder/* 2>/dev/null");
    }
}

#################################################################
######################## ANI subroutines ########################
#################################################################

########## our blast+ reinterpretation
sub calcANIx {
    my($qfile,$sfile) = @_;
    my $qLarge = inflateFile("$qfile");
    my $dbfile = nameDB("$sfile");
    formatDB("$sfile","$dbfile","blast");
    my $action = "plusblasting";
    my $addSpace = " " x ( length($action) - 3 );
    print "$action $qfile\n${addSpace} vs $dbfile\n";
    my @table = qw(
                      6
                      qaccver
                      saccver
                      pident
                      qstart
                      qend
                      length
              );
    #### megablast is the default task, but we're making it explicit
    my $blastrun
        = qq(blastn -query $qLarge -db $dbfile)
        . qq( -outfmt ') . join(" ",@table) . qq(')
        . qq( -evalue $maxEvalue )
        . qq( -dust no )
        . qq( -xdrop_gap 150 )
        . qq( -penalty -1 )
        . qq( -reward 1 )
        . qq( -gapopen 5 )
        . qq( -gapextend 2 )
        . qq( -num_threads $cpus )
        . qq( -task megablast )
        ;
    my $totalLn = 0;
    my $totalId = 0;
    my @acceptedCoords  = ();
  XLINE:
    for my $line ( qx($blastrun) ) {
        chomp($line);
        my($qid,$sid,$pident,
           $qstart,$qend,
           $ln) = split(/\t/,$line);
        next XLINE if( $ln < $cutLn );
        $totalLn += $ln;
        $totalId += ($ln * ($pident/100));
    }
    if( $totalLn > 0 ) {
        my $ani = sprintf("%.${round}f",100*($totalId/$totalLn));
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

########## ANI as described in article (Kostas program), but with blast+
sub calcANIp {
    my($query,$sfile) = @_;
    my $qfile  = namePieces("$query");
    my $dbfile = nameDB("$sfile");
    formatDB("$sfile","$dbfile","blast");
    my $action = "plusblasting";
    my $addSpace = " " x ( length($action) - 3 );
    print "$action $qfile\n${addSpace} vs $dbfile\n";
    my @table = qw(
                      6
                      qaccver
                      saccver
                      pident
              );
    #                  nident
    my $blastrun
        = qq(blastn -query $qfile -db $dbfile)
        . qq( -outfmt ') . join(" ",@table) . qq(')
        . qq( -evalue $maxEvalue )
        . qq( -dust no )
        . qq( -max_hsps 1 )
        . qq( -max_target_seqs 1 )
        . qq( -qcov_hsp_perc 70 )
        . qq( -perc_identity 70 )
        . qq( -xdrop_gap 150 )
        . qq( -penalty -1 )
        . qq( -reward 1 )
        . qq( -gapopen 5 )
        . qq( -gapextend 2 )
        . qq( -num_threads $cpus )
        . qq( -task blastn )
        ;
    #    . qq( -perc_identity 30 )
    my @pcids = ();
    for my $line ( qx($blastrun) ) {
        chomp($line);
        #my($qid,$sid,$pident,$ident) = split(/\s+/,$line);
        my($qid,$sid,$pident) = split(/\t/,$line);
        push(@pcids,$pident);
        ####### if we wanted something as described in paper:
        # my $pcid = 100 * ( $ident / $cutLn );
        # if( $pcid > 30 ) {
        #     push(@pcids,$pcid);
        # }
    }
    my $count = @pcids;
    my $total = sum(@pcids);
    if( $count > 0 ) {
        my $ani = sprintf("%.${round}f",($total/$count));
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

######## Intended to follow article's description
######## but ended filtered as in Kostas' perl program
sub calcANIb {
    my($query,$sfile) = @_;
    my $qfile  = namePieces("$query");
    my $dbfile = nameDB("$sfile");
    formatDB("$sfile","$dbfile","blast");
    my $action = "legacy blasting";
    my $addSpace = " " x ( length($action) - 3 );
    print "$action $qfile\n${addSpace} vs $dbfile\n";
    my $blastrun
        = qq(blastall -p blastn -i $qfile -d $dbfile -e $maxEvalue)
        . qq( -F F )
        . qq( -v 1 )
        . qq( -b 1 )
        . qq( -X 150 )
        . qq( -q -1 )
        . qq( -a $cpus )
        . qq( -m 8 )
        ;
    my $thrLn  = $cutLn * 0.7;
    #print "cutting unless length > $thrLn\n";
    my %pident = ();
    for my $line ( qx($blastrun) ) {
        chomp($line);
        my($qaccver,$saccver,
           $pident,$length,$mismatch,$gapopen,
           $qstart,$qend,
           $sstart,$send,
           $evalue,$bitscore) = split(/\t/,$line);
        ######## filtered as in Kostas' ANI program
        if( $length > $thrLn ## 700 in Kostas' program
                && $pident > 70
                && $pident > $pident{"$qaccver"} ) {
            $pident{"$qaccver"} = $pident;
        }
    }
    my @pidents = values %pident;
    my $total = sum(@pidents);
    my $count = @pidents;
    if( $count > 0 ) {
        my $ani = sprintf("%.${round}f",($total/$count));
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

sub calcANIm {
    my($qfile,$sfile) = @_;
    print "nucmering $qfile vs $sfile\n";
    my $outfile = nakedName("$qfile") . "." . nakedName("$sfile");
    $outfile = join("/",$tempFolder,$outfile);
    my $qLarge = inflateFile("$qfile");
    my $sLarge = inflateFile("$sfile");
    # implemented with "dnadiff" as per:
    # https://doi.org/10.1093/bioinformatics/btv681
    my $nucmerer
        = qq(nucmer --mum -t $cpus -p $outfile $qLarge $sLarge);
    my $outNum  = qx($nucmerer);
    my $dnaDiffer
        = qq(dnadiff -d $outfile.delta -p $outfile);
    my $outDiff = qx($dnaDiffer);
    my $totalLn = 0;
    my $totalId = 0;
    my $ani     = 0;
    my $getANI = 'no';
    open( my $DIFF,"<","$outfile.report" );
  READDIFF:
    while(<$DIFF>) {
        if( $getANI eq 'yes' && m{^AvgIdentity} ) {
            my($label,$ani1,$ani2) = split;
            $ani = sprintf("%.${round}f",( ( $ani1 + $ani2 ) / 2 ) );
            last READDIFF;
        }
        if( m{^M-to-M\s+} ) {
            $getANI = 'yes';
        }
    }
    close($DIFF);
    if( $ani > 0 ) {
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

sub calcANIl {
    my($query,$sfile) = @_;
    my $qfile = namePieces("$query");
    my $dbfile = nameDB("$sfile");
    formatDB("$sfile","$dbfile","lastal");
    my $action = "lastaling";
    my $addSpace = " " x ( length($action) - 3 );
    print "$action $qfile\n${addSpace} vs $dbfile\n";
    my $lastrun
        = qq(lastal)
        . qq( -N 1 )
        . qq( -f BlastTab )
        . qq( -P $cpus )
        . qq( -q 1 )
        . qq( $dbfile $qfile )
        ;
    my $thrLn  = $cutLn * 0.7;
    #print "cutting unless length > $thrLn\n";
    my %pident = ();
  LASTLINE:
    for my $line ( qx($lastrun) ) {
        chomp($line);
        next LASTLINE if( m{^#} );
        my($qaccver,$saccver,
           $pident,$length,$mismatch,$gapopen,
           $qstart,$qend,
           $sstart,$send,
           $evalue,$bitscore) = split(/\t/,$line);
        next LASTLINE if( $evalue > $maxEvalue );
        ######## filtered as in 2007 ANI paper
        if( $length > $thrLn ## 700 in paper
                && $pident > 30
                && $pident > $pident{"$qaccver"} ) {
            $pident{"$qaccver"} = $pident;
        }
    }
    my @pidents = values %pident;
    my $total = sum(@pidents);
    my $count = @pidents;
    if( $count > 0 ) {
        my $ani = sprintf("%.${round}f",($total/$count));
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

sub calcANIu {
    my($query,$sfile) = @_;
    my $qfile = namePieces("$query");
    my $sLarge = inflateFile("$sfile");
    my $action = "ublasting";
    my $addSpace = " " x ( length($action) - 3 );
    print "$action $qfile\n${addSpace} vs $sLarge\n";
    my @table = qw(
                      query
                      target
                      pctpv
              );
    #                  pv
    my $outfile = nakedName("$query") . "." . nakedName("$sfile");
    $outfile = join("/",$tempFolder,$outfile);
    my $ublaster
        = qq(usearch -ublast $qfile -db $sLarge -quiet )
        . qq( -qmask none -dbmask none )
        . qq( -evalue $maxEvalue -accel 0.5 )
        . qq( -query_cov 0.70 )
        . qq( -id 0.30 )
        . qq( -threads $cpus )
        . qq( -strand both )
        . qq( -userout $outfile )
        . qq( -userfields ) . join("+",@table);
    my $outErr = qx($ublaster);
    my %pcid = ();
    open( my $UBLASTED,"<","$outfile" );
    while(<$UBLASTED>) {
        chomp;
        my($qid,$sid,$pcid) = split(/\t/,$_);
        if( $pcid > $pcid{"$qid"} ) {
            $pcid{"$qid"} = $pcid;
        }
    }
    close($UBLASTED);
    my @pcids = values %pcid;
    my $total = sum(@pcids);
    my $count = @pcids;
    if( $count > 0 ) {
        my $ani = sprintf("%.${round}f",($total/$count));
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

sub calcANIf {
    my($qfile,$sfile) = @_;
    print "fastANIing $qfile vs $sfile\n";
    my $outfile = nakedName("$qfile") . "." . nakedName("$sfile");
    $outfile = join("/",$tempFolder,$outfile);
    my $qLarge = inflateFile("$qfile");
    my $sLarge = inflateFile("$sfile");
    my $fastANIer
        = qq(fastANI -o $outfile -q $qLarge -r $sLarge)
        . qq( -t $cpus --fragLen $cutLn )
        ;
    ##### trusting that the fastANI was optimized for length:
    #    . qq( --fragLen $cutLn )
    my $fastOut = qx($fastANIer 2>&1);
    my $anif = 0;
    open( my $FASTANIED,"<","$outfile" );
    while(<$FASTANIED>) {
        chomp;
        my ($gnm1,$gnm2,$ani,@numbers) = split;
        $anif = $ani;
    }
    close($FASTANIED);
    #unlink("$outfile");
    if( $anif > 0 ) {
        my $ani = sprintf("%.${round}f",$anif);
        print "ANI is $ani\n";
        return($ani);
    }
    else {
        return();
    }
}

#################################################################
#################### other subroutines ##########################
#################################################################
sub cutIntoPieces {
    my($inSeqFile,$outSeqFile) = @_;
    if( -f "$outSeqFile" ) {
        print "  the genome is already cut\n";
    }
    else {
        print "  slurping sequences\n";
        my($opener) = how2open("$inSeqFile");
        my %seq_of = ();
        my $id     = '';
        open( my $INFILE,"-|","$opener $inSeqFile" );
        while(<$INFILE>) {
            if( m{>(\S+)} ) {
                $id = $1;
                #print "The ID is ",$id,"\n";
            }
            else {
                chomp;
                $seq_of{"$id"} .= $_;
                #print length($seq_of{"$id"}),"\n";
            }
        }
        close($INFILE);
        print "  cutting sequences\n";
        my $tcount = 0;
        open( my $OUTSEQFILE,">","$outSeqFile" );
        for my $id ( keys %seq_of ) {
            #print "working with ",$id,"\n";
            my $seq   = $seq_of{"$id"};
            my $count = 0;
            my $piece = '';
            while ( length($seq) >= $cutLn ) {
                $seq =~ s{^(\S{$cutLn})}{} and $piece = $1;
                ##### work only with DNA pieces
                if( $piece =~ m{^[ACGT]+$}i ) {
                    print {$OUTSEQFILE} ">",$id,".",$count,"\n",$piece,"\n";
                    $count++;
                }
            }
            $tcount += $count;
        }
        close($OUTSEQFILE);
        print "  printed $tcount ${cutLn}bp long pieces\n";
    }
}

sub nakedName {
    my $inName  = $_[0];
    my $outName = $inName;
    $outName =~ s{^\S+/}{};
    $outName =~ s{\.(Z|gz|bz2)$}{};
    $outName =~ s{\.fna\S*}{};
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

sub namePieces {
    my $fnaFileName = $_[0];
    my $outName = nakedName("$fnaFileName");
    if( length($outName) > 0 ) {
        return("$tempFolder/$outName.pieces");
    }
    else {
        print "there's a problem with $fnaFileName -> $outName\n";
        return();
    }
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
    my( $opener ) = how2open($origfile);
    my $rootName = nakedName("$origfile");
    my $inflated = "$tempFolder/$rootName.fna";
    unless( -f "$inflated" ) {
        system("$opener $origfile > $inflated");
    }
    return("$inflated");
}

sub formatDB {
    my($file,$dbfile,$dbType) = @_;
    my($opener) = how2open("$file");
    if( $dbType eq "blast" ) {
        if( -f "$dbfile.nsq" ) {
            print "  the blast DB file is already there\n";
        }
        else {
            print "producing blastDB: $dbfile\n";
            my $mkblastdb
                = qq($opener $file |)
                . qq( makeblastdb -parse_seqids )
                . qq( -dbtype nucl )
                . qq( -title $dbfile )
                . qq( -out $dbfile );
            system("$mkblastdb 1>/dev/null");
            #system("makembindex -input $dbfile");
        }
    }
    else { ### lastal for now
        if( -f "$dbfile.bck" ) {
            print "  the lastal DB is already there\n";
        }
        else {
            print "producing lastDB: $dbfile\n";
            my $mklastdb
                = qq( $opener $file | )
                . qq(lastdb $dbfile );
            system("$mklastdb 1>/dev/null");
        }
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

sub findFnaFiles {
    my($testDir,$label) = @_;
    print "finding $label fna|fasta files in $testDir\n";
    my @foundFiles = ();
    #### now repopulate arrays:
    opendir( my $FNADIR,"$testDir" );
    my @fnaFiles = sort grep { m{^\w} && m{\.(fna|fasta)} } readdir $FNADIR;
    for my $fnaFile ( @fnaFiles ) {
        my $cntSeqs = checkFasta("$testDir/$fnaFile");
        if( $cntSeqs > 0 ) {
            push( @foundFiles,"$testDir/$fnaFile" );
        }
        else {
            print " ignoring empty file: $fnaFile\n";
        }
    }
    my $cntFound = @foundFiles;
    if( $cntFound > 0 ) {
        print " will work with $cntFound $label files\n";
        return(@foundFiles);
    }
    else {
        die " no fna/fasta files in $label directory ($testDir)\n\n";
    }
}

sub findWorkFiles {
    my ($item,@testFiles) = @_;
    my $toFind = @testFiles;
    print "checking for $toFind $item files\n";
    my @found = ();
    for my $file ( sort @testFiles ) {
        if( -f "$file" ) {
            my $cntSeqs = checkFasta("$file");
            if( $cntSeqs > 0 ) {
                push(@found,$file);
            }
            else {
                print " $file does not have sequences\n";
            }
        }
        else {
            print " did not find $file\n";
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

sub checkFasta {
    my $checkFasta = $_[0];
    my( $opener,$file ) = how2open("$checkFasta");
    my $cntSeqs = 0;
    open( my $CFA,"-|","$opener $file" );
    while(<$CFA>) {
        if( m{^>} ) {
            $cntSeqs++;
        }
    }
    close($CFA);
    return($cntSeqs);
}

