#!/usr/bin/perl

#########################################################
#                                                       #
#       Author : Gabo Moreno-Hagelsieb                  #
#                                                       #
#########################################################

use strict;
use Getopt::Long;
use Pod::Text;
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);
my $columns = qx(tput cols);
chomp($columns);
my $width = $columns >= 100 ? 80 : $columns - 3;
my $parser = Pod::Text->new (sentence => 0, width => $width, margin => 1);

my $ownName = $0;
$ownName =~ s{.*/}{};

my @famDBs = qw(
                   cog
                   cdd
                   pfam
                   tigrfam
                   superfamily
                   VFDB
                   Toxins
           );
my $matchDBs = join('|',@famDBs);

my @programs = qw(
                     hmmscan
                     rpsblast
                     mmseqs
             );

my %defaultProg = (
    'cog'         => 'rpsblast',
    'cdd'         => 'rpsblast',
    'pfam'        => 'hmmscan',
    'tigrfam'     => 'hmmscan',
    'superfamily' => 'hmmscan',
    'vfdb'        => 'hmmscan',
    'toxins'      => 'hmmscan',
);

my $defCPU      = 2;
my $matchProg   = join('|',@programs);
my $defEvalue   = 1e-3;
my $ncbiEvalue  = 1e-2;
my $defmmAltAli = 10;
my @queries     = ();
my $famDB       = '';
my $scanProgram = '';
my $resultsDir  = 'scanFams';
my $cpus        = $defCPU;
my $Evalue      = '';
my $overwrite   = 'F';
my $cluster     = 'F';
my $baseDBdir
    = -d $ENV{"FAMDB"}                ? $ENV{"FAMDB"}
    : -d "/usr/local/famDBs"          ? "/usr/local/famDBs"
    : -d "/usr/local/domainDBs"       ? "/usr/local/domainDBs"
    : -d "/usr/local/DB/domainDBs"    ? "/usr/local/DB/domainDBs"
    : -d "/ResearchData/DB/domainDBs" ? "/ResearchData/DB/domainDBs"
    : -d $ENV{"GENOMEDB"}             ? $ENV{"GENOMEDB"} . "/domainDBs"
    : ".";
### choosing database
my $defaultNA = "none";
my $cddDB
    = -d "$baseDBdir/cddDB" ? "$baseDBdir/cddDB"
    : -d "cddDB"            ? "cddDB"
    : $defaultNA;
my $xfamDB
    = -d "$baseDBdir/xfamDB" ? "$baseDBdir/xfamDB"
    : -d "xfamDB"            ? "xfamDB"
    : $defaultNA;
my $mmseqsDB
    = -d "$baseDBdir/mmseqsDB" ? "$baseDBdir/mmseqsDB"
    : -d "mmseqsDB"            ? "mmseqsDB"
    : $defaultNA;

my %progDBdir = (
    'rpsblast' => "$cddDB",
    'hmmscan'  => "$xfamDB",
    'mmseqs'   => "$mmseqsDB",
);

### check if there's more than one processor or assume there's 2.
my $cpu_count
    = qx(getconf _NPROCESSORS_ONLN 2>/dev/null)
    =~ m{\s*(\d+)\s*}                 ? $1
    : qx(sysctl -a 2>/dev/null | grep 'cpu.thread_count')
    =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
    =~ m{\.max-threads\s+=\s+(\d+)}    ? $1
    : $defCPU;

my $podUsage
    = qq(=pod\n\n)
    . qq(=head1 NAME\n\n)
    . qq($ownName - Scan protein sequences against protein family profiles\n\n)
    . qq(=head1 SYNOPSIS\n\n)
    . qq($ownName -q [fastaFile] -f [$matchDBs] [options]\n\n)
    . qq(=head1 EXAMPLE\n\n)
    . qq($ownName -q GCF_000005845.faa.gz -f pfam -p mmseqs -o PFAM\n\n)
    . qq($ownName -q fastaFiles/*.faa.gz -f cdd -o CDD\n\n)
    . qq(=head1 OPTIONS\n\n)
    . qq(=over\n\n)
    . qq(=item B<-q>\n\n)
    . qq(query fasta file(s), required\n\n)
    . qq(=item B<-f>\n\n)
    . qq(family database: file or [$matchDBs], required\n\n)
    . qq(=item B<-p>\n\n)
    . qq(program [$matchProg], except for mmseqs, can be guessed from\n)
    . qq(database:\n\n)
    . qq(=over\n\n)
    . qq(=item -\n\n)
    . qq(rpsblast for cog and cdd\n\n)
    . qq(=item -\n\n)
    . qq(hmmscan for Pfam and TIGRFAM\n\n)
    . qq(=item -\n\n)
    . qq(mmseqs has to be specified in command line\n\n)
    . qq(=back\n\n)
    . qq(=item B<-e>\n\n)
    . qq(e-value threshold, default $defEvalue (NCBI uses $ncbiEvalue),\n)
    . qq(scientific notation acceptable (e.g. 1e-3)\n\n)
    . qq(=item B<-o>\n\n)
    . qq(output folder, default: $resultsDir\n\n)
    . qq(=item B<-w>\n\n)
    . qq(overwrite existing result files [T|F]: default $overwrite\n\n)
    . qq(=item B<-x>\n\n)
    . qq(number of CPUs to use, default: 2 (max: $cpu_count)\n\n)
    . qq(=item B<-c>\n\n)
    . qq(running in computer cluster [T|F], default 'F'\n\n)
    . qq(=back\n\n)
    . qq(=head1 DESCRIPTION\n\n)
    . qq(This program scans protein sequences against [$matchDBs]\n)
    . qq( databases using either of hmmscan, rpsblast, or mmseqs\n\n)
    . qq(=cut\n\n)
    ;

GetOptions(
    "q=s{,}" => \@queries,
    "f=s"    => \$famDB,
    "p=s"    => \$scanProgram,
    "o=s"    => \$resultsDir,
    "e=f"    => \$Evalue,
    "w=s"    => \$overwrite,
    "x=i"    => \$cpus,
    "c=s"    => \$cluster,
) or podhelp();

if ( !$queries[0] || !$famDB ) {
    podhelp("I need a protein fasta files and a family database:");
}

my @queries = do { my %seen; grep { !$seen{$_}++ } @queries };
my $countQueries = @queries;
my $foundQueries = 0;
my @missing      = ();
my %faaFile      = ();
for my $query ( @queries ) {
    if( -f "$query" ) {
        my $faaFile = $query;
        $query =~ s{\S+/}{};
        $query =~ s{\.(faa|fasta)\S*}{};
        $faaFile{"$query"} = $faaFile;
        $foundQueries++;
    }
    else {
        push(@missing,$query);
    }
}
if( $foundQueries == $countQueries ) {
    print "found $countQueries query fasta files\n";
}
else {
    my $missing = @missing;
    my $error
        = "missing $missing query file(s):\n" . join("\n",@missing) . "\n";
    die $error;
}

$famDB = lc($famDB);
$famDB =~ s{\.(psq|hmm)\S+}{};
$famDB =~ s{^\S+/}{};
$famDB =~ s{\.\S+}{};

my $cpus = $cpus > 0 && $cpus <= $cpu_count ? $cpus : $defCPU;
print "using $cpus cpu threads\n";

### overwrite existing results file?
my $overwrite = $overwrite =~ m{(T|F)}i ? uc($1) : 'F';
print "overwriting results: $overwrite\n";

### working in cluster?
my $cluster = $cluster =~ m{(T|F)}i ? uc($1) : 'F';
print "working in cluster: $cluster\n";

### choosing database
### (will need better ways fo figure out the available databases
###  and avoit genomeTools)

my $scanProgram
    = $scanProgram =~ m{^($matchProg)$}i ? $1
    : exists $defaultProg{"$famDB"}   ? $defaultProg{"$famDB"}
    : "none";

if( $scanProgram eq "none" ) {
    podhelp("I need a program to make these comparisons [$matchProg]")
}

my $Evalue
    = $Evalue > 0 && $Evalue < 1 ? $Evalue
    : $scanProgram eq 'hmmscan' && $famDB =~ m{pfam}i
    ? '--cut_ga'
    : $defEvalue;
print "using an E-value of $Evalue\n";

my $dbPath = $progDBdir{"$scanProgram"};
####### test for family database at the appropriate path:
opendir( my $DBDIR,"$dbPath" );
my @dbFiles
    = sort { -s $b <=> -s $a } grep { m{($famDB)\.}i } readdir($DBDIR);
closedir($DBDIR);
my $dbName = $dbFiles[0] =~ m{($famDB)}i ? $1 : 'none';
my $fullDB = $dbPath . "/" . $dbName;
if( $dbName eq 'none' ) {
    die "there is no $famDB database for $scanProgram at:\n$dbPath\n\n";
}
else {
    if( $scanProgram eq 'hmmscan' ) {
        $fullDB .= ".hmm";
    }
    print "working with database:\n   $fullDB\n";
}

### fields for rpsblast
my @tableFields = qw(
                        qaccver
                        saccver
                        evalue
                        bitscore
                        qstart
                        qend
                        qlen
                        sstart
                        send
                        slen
                 );
#                        qseqid
#                        sseqid

my $tableFields = qq(') . join(" ","6",@tableFields) . qq(');

### fields for mmseqs
my @mmseqsfields = qw(
                         query
                         target
                         evalue
                         bits
                         qstart
                         qend
                         qcov
                         tstart
                         tend
                         tcov
                 );

my $mmseqsfields = join(",",@mmseqsfields);

system "mkdir -p $resultsDir" unless( -d "$resultsDir" );
my $tmpDir = tempdir("/tmp/$ownName.XXXXXXXXXXXX");

my $refName = readNames("$fullDB") if( $scanProgram eq "rpsblast" );

SCANRUNS:
for my $query ( sort keys %faaFile ) {
    my $faaFile = $faaFile{"$query"};
    my ($cat,$queryFile) = figureCompression("$faaFile");
    my $file_name  = join(".",$query,lc($famDB),$scanProgram);
    my $tmpFile    = "$tmpDir/$file_name";
    my $outFile    = "$resultsDir/$file_name";
    my $tmpQuery   = $queryFile;
    $tmpQuery      =~ s{\S+/}{};
    $tmpQuery      = "$tmpDir/$tmpQuery";
    if( -f "$outFile.bz2" ) {
        if( $overwrite eq 'T' ) {
            print "   will overwrite current $outFile.bz2\n";
        }
        else {
            print " won't compare $query vs $famDB\n";
            print "    there's already an $outFile.bz2\n";
            next SCANRUNS;
        }
    }
    print " running $scanProgram $query vs $famDB\n";
    if( $cluster eq 'T' ) {
        print "    preparing database\n";
        system("cp $fullDB* $tmpDir/");
        print "    preparing query\n";
        system("cp $queryFile $tmpQuery");
    }
    my $dbFile
        = $cluster eq 'T' ? "$tmpDir/$famDB"
        : "$fullDB";
    if( $scanProgram eq "rpsblast" ) {
        my $input
            = $cluster eq 'T' ? qq($cat $tmpQuery)
            : qq($cat $queryFile);
        my $rpsblast_command
            = qq( rpsblast -query - )
            . qq( -db $dbFile )
            . qq( -seg yes -soft_masking true )
            . qq( -num_threads $cpus )
            . qq( -evalue $Evalue )
            . qq( -parse_deflines )
            . qq( -comp_based_stats 0 )
            . qq( -outfmt $tableFields );
        my $logStuff
            = qx($input | $rpsblast_command | bzip2 -9 > $tmpFile.bz2 2>&1);
        if( verifyResults("$tmpFile.bz2") ) {
            system("mv $tmpFile.bz2 $outFile.bz2 2>/dev/null");
        }
        else {
            print "   no $famDB matches for $query\n";
        }
    }
    elsif( $scanProgram eq "mmseqs" ) {
        unless( -f "$tmpQuery" ) {
            print "    preparing query\n";
            system("cp $queryFile $tmpQuery");
        }
        my $input = $tmpQuery;
        my $mmseqs_command
            = qq( mmseqs easy-search $input $dbFile $tmpFile $tmpDir)
            . qq( -e $Evalue )
            . qq( --alt-ali $defmmAltAli )
            . qq( --threads $cpus )
            . qq( --comp-bias-corr 0 )
            . qq( --format-output "$mmseqsfields" );
        print "    running mmseqs\n";
        my $logStuff = qx($mmseqs_command 2>&1);
        if( verifyResults("$tmpFile") ) {
            system("mv $tmpFile.bz2 $outFile.bz2 2>/dev/null");
        }
        else {
            print "   no $famDB matches for $query\n";
            my $logFile = "$resultsDir/$scanProgram.log";
            open( my $LOG,">","$logFile" );
            print {$LOG} $logStuff;
            close($LOG);
            print "   output from $scanProgram in $logFile\n";
        }
    }
    else { #### default is hmmscan
        my $input
            = $cluster eq 'T' ? qq($cat $tmpQuery)
            : qq($cat $queryFile);
        my $threshold = $Evalue =~ m{cut_ga} ? $Evalue : "-E $Evalue";
        my $hmmscan_command
            = qq( hmmscan --cpu $cpus --noali $threshold -o /dev/null)
            . qq( --domtblout $tmpFile $dbFile - );
        my $logStuff = qx($input | $hmmscan_command 2>&1);
        if( verifyResults("$tmpFile") ) {
            system("mv $tmpFile*.bz2 $resultsDir/ 2>/dev/null");
        }
        else {
            print "   no $famDB matches for $query\n";
            my $logFile = "$resultsDir/$scanProgram.log";
            open( my $LOG,">","$logFile" );
            print {$LOG} $logStuff;
            close($LOG);
            print "   output from $scanProgram in $logFile\n";
        }
    }
}

if( -d "$tmpDir" ) {
    print "\tcleaning up ...";
    system "rm -rf $tmpDir";
}
print "\n    done with $ownName\n\n";

sub signalHandler {
    if( -d "$tmpDir" ) {
        print "\n\tcleaning up ...";
        system "rm -rf $tmpDir";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n\n";
        die  "    done!\n\n";
    }
}

sub verifyResults {
    my $tmpFile = $_[0];
    my $openTest
        = $tmpFile =~ m{\.bz2$} ? "bzip2 -qdc $tmpFile" : "cat $tmpFile";
    my $verify = 0;
    open( my $TESTF,"-|","$openTest" );
    open( my $VERIFIED,"|-","bzip2 > $tmpFile.tmp" );
    print {$VERIFIED} join("\t",
                           "#Query",
                           "Family",
                           "Evalue",
                           "Bits",
                           "Qstart",
                           "Qend",
                           "Qcov",
                           "Fstart",
                           "Fend",
                           "Fcov"
                       ),"\n";
    while(<$TESTF>) {
        if( $_ !~ m{^#} ) {
            if( $scanProgram eq "mmseqs" ) {
                print {$VERIFIED} $_;
            }
            elsif( $scanProgram eq "rpsblast" ) {
                chomp;
                my( $qseqid,
                    $sseqid,
                    $evalue,
                    $bitscore,
                    $qstart,
                    $qend,
                    $qlen,
                    $sstart,
                    $send,
                    $slen
                ) = split(/\t/,$_);
                $qseqid =~ s{(lcl|ref)\|}{};
                $qseqid =~ s{\|$}{};
                $sseqid =~ s{^(CDD\:|gnl\|CDD\|)}{};
                $sseqid
                    = exists $refName->{"$sseqid"} ? $refName->{"$sseqid"}
                    : $sseqid;
                my $qcov = calcCoverage($qstart,$qend,$qlen);
                my $scov = calcCoverage($sstart,$send,$slen);
                print {$VERIFIED} join("\t",
                                       $qseqid,
                                       $sseqid,
                                       $evalue,
                                       $bitscore,
                                       $qstart,$qend,$qcov,
                                       $sstart,$send,$scov
                                   ),"\n";
            }
            else { ### here default is hmmscan
                my( $dom_name,$dom_id,$dom_len,$qseqid,$caca1,$q_ln,
                    $tevalue,$tscore,$tbias,$n1,$n2,
                    $c_eval,$i_eval,$score,$bias,
                    $domStart,$domEnd,
                    $alnStart,$alnEnd,
                    $qStart,$qEnd,$acc
                ) = split(/\s+/,$_);
                $qseqid =~ s{(lcl|ref)\|}{};
                $qseqid =~ s{\|$}{};
                my $dom
                    = $dom_id eq "-" ? $dom_name
                    : $dom_id;
                my $qcov = calcCoverage($qStart,$qEnd,$q_ln);
                my $dcov = calcCoverage($domStart,$domEnd,$dom_len);
                print {$VERIFIED} join("\t",
                                       $qseqid,
                                       $dom,
                                       $c_eval,
                                       $score,
                                       $qStart,$qEnd,$qcov,
                                       $domStart,$domEnd,$dcov
                                   ),"\n";
            }
            $verify++;
        }
    }
    close($TESTF);
    close($VERIFIED);
    if( $verify > 0 ) {
        if( $scanProgram eq "hmmscan" ) {
            my $newname = $tmpFile;
            $newname =~ s{hmmscan}{hmmscan.original};
            rename("$tmpFile","$newname");
            system qq(bzip2 --best $newname);
            rename("$tmpFile.tmp","$tmpFile.bz2");
        }
        else {
            if( $tmpFile =~ m{\.bz2$} ) {
                rename("$tmpFile.tmp","$tmpFile");
            }
            else {
                rename("$tmpFile.tmp","$tmpFile.bz2");
            }
        }
        return($verify);
    }
    else {
        return();
    }
}

sub podhelp {
    my $extraMessage = $_[0];
    open( my $PIPE,"|-","cat" );
    if( length "$extraMessage" > 2 ) {
        print {$PIPE} "    ",$extraMessage,"\n\n";
    }
    $parser->output_fh($PIPE);
    $parser->parse_string_document($podUsage);
    exit;
}

sub figureCompression {
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

sub calcCoverage {
    my($start,$end,$ln) = @_;
    if( $ln < 1 ) {
        return();
    }
    else {
        my $coverage = ( $end - $start + 1 ) / $ln;
        my $rounded  = sprintf( "%.3f", $coverage );
        return($rounded);
    }
}

sub readNames {
    my $db = $_[0];
    print "learning equiv names for rpsblast:\n $db\n";
    my $dbcmd = qq(blastdbcmd -db $db -entry all -outfmt "%a %t");
    #print $dbcmd,"\n";
    my %name = ();
    for my $line ( qx($dbcmd) ) {
        my($acc,$name,@etc) = split(/\s+/,$line);
        $name =~ s{,$}{};
        if( $name =~ m{^\w+\d+} ) {
            $acc =~ s{^CDD\:}{};
            $name{"$acc"} = $name;
        }
    }
    return(\%name);
}
