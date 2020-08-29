#!/usr/bin/perl

## program to get Ka/Ks (dN/dS) stuff using PAML programs
use strict;
use Getopt::Long;
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my @methods = qw(
                    bayes
                    codeml
                    yn00
            );
my $methods = join("|",@methods);
my @pwProgs = qw(
                    blastp
                    diamond
                    mmseqs
            );
my $matchProg = join("|",@pwProgs);

my $defaultGD = 'GenomeData';
my $maxEvalue = 1e-6;
my $defPW     = 'diamond';

############ variables:
my $queryGnm   = '';
my $targetGnm = '';
my $genomeData = $defaultGD;
my $method     = '';
my $codonF     = 2;
my $outDir     = '';
my $pwProg     = $defPW;
my $pwDir      = "compRuns";

my $ownName = $0;
$ownName =~ s{.*/}{};

my $helpBrief
    = qq(about:\n)
    . qq(    This program produces files with dN/dS results\n\n)
    . qq(usage:\n)
    . qq(    $ownName -q <queryFile> -t <targetFile> [options]\n)
    . qq(\noptions:\n)
    . qq(   -q query genome ID [GCF_000005845], required\n)
    . qq(   -t target genome ID [GCF_000287335], required\n)
    . qq(   -g directory with genome data (faa, cds|ffn),\n)
    . qq(       default $defaultGD\n)
    . qq(   -w dN/dS estimate method [$methods], default bayes\n)
    . qq(   -c codon freq estimate, 2: Goldman Yang, 5: Muse Gaut,\n)
    . qq(      default 2 (Goldman Yang)\n)
    . qq(   -p program for pairwise comparisons [$matchProg],\n)
    . qq(       default: $defPW\n)
    . qq(       (CAUTION: mmseqs is not working yet)\n)
    . qq(   -m directory for $matchProg results, default compRuns\n)
    . qq(   -o output directory, default [method]-dNdS\n)
    . qq(\n)
    . qq(requirements:\n)
    . qq(    This program requires the blastp program, or blast results\n)
    . qq(    for the query + target genomes in the compRuns directory,\n)
    . qq(    or in the directory indicated with the '-m' option\n\n)
    . qq(    The program also requires at least two files per genome:\n)
    . qq(       1. A fasta file with protein sequences ending in ".faa"\n)
    . qq(       2. A fasta file with corresponding coding sequences\n)
    . qq(          ending in ".cds" or ".ffn"\n)
    . qq(    The files can be compressed with gzip or bzip2\n\n)
    ;

GetOptions(
    "q=s" => \$queryGnm,
    "t=s" => \$targetGnm,
    "g=s" => \$genomeData,
    "w=s" => \$method,
    "c=s" => \$codonF,
    "o=s" => \$outDir,
    "p=s" => \$pwProg,
    "m=s" => \$pwDir,
) or die $helpBrief;


if( !$queryGnm || !$targetGnm ) {
    print $helpBrief;
    exit;
}

my $method  = $method =~ m{^($methods)$}i    ? lc($method) : "bayes";
my $program = $method =~ m{^(bayes|codeml)$} ? "codeml"    : "yn00";
my $codonF  = $codonF =~ m{^(2|5)$}          ? $codonF     : 2;
my $testData = ( -d "$genomeData" )          ? $genomeData : "NA";
if( $testData eq "NA" ) {
    die "There's no genome data directory ($genomeData)\n\n".$helpBrief;
}
else {
    print "Genome data to be found at $genomeData\n";
}

#### directories where to find files:
my $cwd = qx(pwd);
chomp($cwd);
my $workDir = tempdir("/tmp/dNdS.XXXXXXXXXXXX");
my $outDir  = length("$outDir") > 1 ? $outDir : "$method-dNdS";
unless( -d "$outDir" ) {
    mkdir("$outDir");
}
unless( -d "$outDir/$queryGnm" ) {
    mkdir("$outDir/$queryGnm");
}
my $outFile  = $outDir  . "/$queryGnm/$targetGnm.$method";
my $tmpFile  = $workDir . "/$queryGnm.$targetGnm.$method";
my $bugFile  = $outDir  . "/$queryGnm/$targetGnm.$method.problem";
my $bugTmp   = $workDir . "/$queryGnm.$targetGnm.problem";
my $dNdSres  = "$queryGnm.$targetGnm.dNdS";
my $alnseq   = "$queryGnm.$targetGnm.nuc";
####### in case we're running pairwise comparisons:
my $pwProg   = $pwProg =~ m{^($matchProg)$} ? lc($1) : $defPW;
my $pwDB     = $workDir . "/$targetGnm";
my $alnFile  = $pwDir   . "/$queryGnm/$targetGnm." . $pwProg . ".bz2";

#### check if we have pairwise comparison results
#### run pw comparison otherwise
if( -f "$alnFile" ) {
    unless( my $good2go = checkAlignSeqs($alnFile) ) {
        print "comparing sequences with $pwProg\nwill save to $alnFile\n";
        runPW("$genomeData/$queryGnm.faa",
              "$genomeData/$targetGnm.faa",
              "$pwProg");
    }
}
else {
    print "comparing sequences with $pwProg\nwill save to $alnFile\n";
    runPW("$genomeData/$queryGnm.faa",
          "$genomeData/$targetGnm.faa",
          "$pwProg");
}

#### control files
### yn00 and codeml need a special control file.
### This file determines how the program runs and how it calculates
### the dN, dS and dN/dS values.
produceControlFile("$program");

## then we need to read the dna sequences coding each protein sequence,
## cut in triplets (codons), and relate each codon to its aminoacid
my $refSeq = readFastaSeqs("$genomeData","$queryGnm","$targetGnm") or
    die "    some CDS or FFN file was unavailable\n\n";

## then we need to read the alignment sequences
## and calculate the dN/dS values
## this is the result of blastp with a tab-separated
## table output format

### read blast results and check the number of pairs to run through dNdS
### analyses
my ($countLines,$refAlignments) = readPW($alnFile);

### now calculate dN/dS
if( $countLines > 0 ) {
    print "calculating dNdS for $countLines pairs\n";
    my $countDone = 0;
    my $countBugs = 0;
    open( our $DNDS,">","$tmpFile");
    if( $method eq "bayes" ) {
        print {$DNDS} join("\t","Gene_i","Gene_j",
                           "t","t(se)","t(ml)","N","S",
                           "dN","dS","dN/dS","dN/dS(se)","dN/dS(ml)","P1"),"\n";
    }
    elsif( $method eq "codeml" ) {
        print {$DNDS} join("\t","Gene_i","Gene_j",
                           "t","N","S",
                           "dN","dS","dN/dS"),"\n";
    }
    else { ## $method eq "yn00"
        print {$DNDS} join("\t","Gene_i","Gene_j",
                           "t","N","S",
                           "dN","dN(se)","dS","dS(se)","dN/dS"),"\n";
    }
    ### file for problems
    open( our $BUGS,">","$bugTmp");
    ### setup feedback counter
    my $checkpnt = 0;
    ### end setup feedback counter
    for my $alignment ( @{$refAlignments} ) {
        my (
            $query,$target,
            $qstart,$qend,
            $sstart,$send,
            $qalnseq,$salnseq
        ) = split(/\t/,$alignment);
        if( my $DNAstatus
                = produceDNAalignment(
                    $query,$target,
                    $qstart,$qend,
                    $sstart,$send,
                    $qalnseq,$salnseq
                ) ) {
            if( $DNAstatus eq "equal" ) {
                print {$BUGS} join("\t",$query,$target,
                                   "identical DNA sequences"),"\n";
                $countBugs++;
            }
            else {
                my ($tt,$tN,$tS,$dN,$dS,$dNdS) = getdNdS();
                #print join(", ",$tt,$tN,$tS,$dN,$dS,$dNdS),"<-after\n";exit;
                #print join("\t",$query,$target,$dn,$ds,$dnds),"\n";
                ##### making sure that we've got a value from this thing
                if( $dNdS =~ m{^\d+} ) {
                    my @tt = split(/\t/,$tt);
                    my @dS = split(/\t/,$dS);
                    my $problem = "";
                    ####### Bayes seems to always work
                    if( $method ne "bayes" ) {
                        if( $tt[0] > 20 ) {
                            $problem .= "saturation";
                        }
                        elsif( $dS[0] <= 0 ) {
                            $problem .= "unchanged";
                        }
                    }
                    if( length("$problem") > 2 ) {
                        print {$BUGS}
                            join("\t",$query,$target,$problem.": "),
                            join(", ",$tt,$tN,$tS,$dN,$dS,$dNdS),"\n";
                        $countBugs++;
                    }
                    else {
                        print {$DNDS}
                            join("\t",
                                 $query,$target,
                                 $tt,$tN,$tS,$dN,$dS,$dNdS
                             ),"\n";
                        $countDone++;
                    }
                }
                else {
                    print {$BUGS}
                        join("\t",$query,$target,"unknown: "),
                        join(", ",$tt,$tN,$tS,$dN,$dS,$dNdS),"\n";
                    $countBugs++;
                }
            }
        }
        else {
            print {$BUGS}
                join("\t",$query,$target,"failed DNA alignment"),"\n";
            $countBugs++;
        }
        ### feedback to terminal
        $checkpnt++;
        progressLine($checkpnt,$countLines,4);
        ### end feedback
    }
    close($DNDS);
    close($BUGS);
    if( $countBugs > 0 ) {
        system("bzip2 -9 $bugTmp");
        system qq(mv $bugTmp.bz2 $bugFile.bz2 &>/dev/null);
    }
    if( $countDone > 0 ) {
        print "\tobtained $countDone dN/dS values\n";
        system("bzip2 -9 $tmpFile");
        system qq(mv $tmpFile.bz2 $outFile.bz2 &>/dev/null);
    }
    else {
        if( $countBugs > 0 ) {
            system("bzip2 -9 $bugTmp");
            system qq(mv $bugTmp.bz2 $bugFile.bz2 &>/dev/null);
            print "\n\tseems like we had no results. Check:\n"
                . "  $bugFile.bz2\n";
        }
        else {
            print "\n\tseems like we had no results\n"
                . "  no $bugFile produced either\n";
        }
    }
}
else {
    print "no alignments suitable for dN/dS calculations\n";
}

##### clean up
if( -d "$workDir" ) {
    print "\tcleaning up ...\n";
    system "rm -r $workDir";
}
print "\tdone!\n\n";

#################################################################
######################## subroutines ############################
#################################################################
sub getdNdS {
    my $tt   = "-";
    my $tN   = "-";
    my $tS   = "-";
    my $dS   = "-";
    my $dN   = "-";
    my $dNdS = "-";
    my $p1   = "-"; # to save prob > 1
    ## we need to run the program in the same place where the control
    ## and the alignments files are
    system <<"RUNMETHOD";
cd $workDir
$program <<FIN >& /dev/null

FIN

cd $cwd

RUNMETHOD
    ### now there should be a results file
    if( -s "$workDir/$dNdSres" ) {
        #print "opening $workDir/$dNdSres ($method)\n";exit;
        open( my $KK,"<","$workDir/$dNdSres");
        my @lines = (<$KK>);
        close($KK);
        ##### first program codeml (method can be codeml or bayes):
        if( $program eq "codeml" ) {
            if( my @matching = grep { m{^t=\s+\d\S+\s+S\s*=} } @lines ) {
                my $test = $matching[-1];
                $tt   = $test =~ m{^t\s*=\s*(\d\S+)}      ? $1 : "-";
                $tN   = $test =~ m{\sN\s*=\s*(\d\S+)}     ? $1 : "-";
                $tS   = $test =~ m{\sS\s*=\s*(\d\S+)}     ? $1 : "-";
                $dN   = $test =~ m{\sdN\s*=\s*(\d\S+)}    ? $1 : "-";
                $dS   = $test =~ m{\sdS\s*=\s*(\d\S+)}    ? $1 : "-";
                $dNdS = $test =~ m{\sdN/dS\s*=\s*(\d\S+)} ? $1 : "-";
            }
            if( $method eq "bayes" ) {
                ### the last line contains the Bayesian calculations
                my $lastLine = $lines[-1];
                $lastLine =~ s{^\s+}{};
                if( $lastLine =~ m{^\d+} ) {
                    my ( $bT,$bW,$seT,$seW,$covTW,$corrTW,$pW1 )
                        = split(/\s+/,$lastLine);
                    $tt   = join("\t",$bT,$seT,$tt);
                    $dNdS = join("\t",$bW,$seW,$dNdS,$pW1);
                }
                else {
                    $tt   = join("\t","NA","NA",$tt);
                    $dNdS = join("\t","NA","NA",$dNdS,"NA");
                }
            }
        }
        elsif( $program eq "yn00" ) {
            if( my @matching = grep { m{\d+\s*\+\-} } @lines ) {
                my $test = $matching[-1];
                $test =~ s{^\s+}{};
                my @a = split(/\s+/,$test);
                $tt   = $a[4];
                $tN   = $a[3];
                $tS   = $a[2];
                $dN   = join("\t",$a[7],$a[9]);
                $dS   = join("\t",$a[10],$a[12]);
                $dNdS = $a[6];
            }
        }
    }
    if( -f "$workDir/$dNdSres" ){
        unlink("$workDir/$dNdSres");
    }
    unlink("$workDir/$alnseq");
    #print join(", ",$tt,$tN,$tS,$dN,$dS,$dNdS),"<-before\n";
    return($tt,$tN,$tS,$dN,$dS,$dNdS);
}

sub produceDNAalignment {
    my( $query,$target,
        $qstart,$qend,
        $sstart,$send,
        $qalnseq,$salnseq
    ) = @_;
    my @qaln = split(//,$qalnseq);
    my @saln = split(//,$salnseq);
    ## get codons for query alignment
    my @qcodons = getCodons($refSeq->{"$query"});
    ## get codons for target alignment
    my @scodons = getCodons($refSeq->{"$target"});
    my $qNcodons = @qcodons;
    my $sNcodons = @scodons;
    if( $qNcodons < 10 ) {
        return undef;
    }
    if( $sNcodons < 10 ) {
        return undef;
    }
    else {
        ## now it is time to produce the dna aligned by
        ## codons following the protein alignment
        my $qoffset = $qstart - 1;
        my $soffset = $sstart - 1;
        my $qalign  = "";
        my $salign  = "";
        for my $n ( 0 .. $#qaln ) {
            if( $qaln[$n] =~ /\w/ ) {
                $qalign .= $qcodons[$qoffset];
                if( $saln[$n] =~ /\w/ ) {
                    $salign .= $scodons[$soffset];
                    $soffset++;
                }
                else {
                    $salign .= "---";
                }
                $qoffset++;
            }
            else {
                $qalign .= "---";
                if($saln[$n] =~ /\w/) {
                    $salign .= $scodons[$soffset];
                    $soffset++;
                }
                else {
                    ## nothing happens for this one
                }
            }
        }
        if( "$qalign" eq "$salign" ) {
            return("equal");
        }
        elsif(
            length($qalign) > 0 &&
                length($salign) > 0 &&
                length($qalign) eq length($salign)
            ) {
            open( my $FNAK,">","$workDir/$alnseq" );
            print {$FNAK}
                "\t2\t",
                join("\n",length($qalign),$query,$qalign,$target,$salign),"\n";
            close($FNAK);
            #print "alignment printed:\n   $workDir/$alnseq\n";exit;
            return("printed");
        }
        else {
            return undef;
        }
    }
}

sub signalHandler {
    if( -d "$workDir" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $workDir";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n";
        die  "    done!\n\n";
    }
}

sub produceControlFile {
    my $program = $_[0];
    my $controlFile = "$workDir/$program.ctl";
    my $runmode = $method eq "bayes" ? "-3" : "-2";
    if( $program eq "codeml" ) {
        open( my $CTLC,">","$controlFile");
        print {$CTLC} <<"FINC";
      seqfile = $alnseq * sequence data file name
      outfile = $dNdSres * main result file name
        noisy = 0      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = $runmode     * -2:pairwise; -3 pairwise Bayesian
      seqtype = 1      * 1:codons
    CodonFreq = $codonF      * 0:equal, 1:F1X4, 2:F3X4, 3:F61, 5:Bielawski
        model = 0      *
      NSsites = 0      * 0:one, 2:selection
        icode = 0      * 0:universal code
    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * fixed or initial value
    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated
        omega = 0.4    * initial omega value
FINC
        close($CTLC);
    }
    elsif( $program eq "yn00" ) {
        open( my $CTLY,">","$controlFile");
        print {$CTLY} <<"FINY";
      seqfile = $alnseq      * sequence data file name
      outfile = $dNdSres     * main result file
      verbose = 0 * 1: detailed output (list sequences), 0: concise output
        icode = 0 * 0:universal code; 1:mammalian mt; 2-10:see below
    weighting = 0 * weighting pathways between codons (0/1)?
   commonf3x4 = 0 * use one set of codon freqs for all pairs (0/1)?
FINY
        close($CTLY);
    }
    print "control file written:\n   $controlFile\n";
}

sub readFastaSeqs {
    my($fastaDir,@gnms) = @_;
    my %seq = ();
    my $cdsMistakes = 0;
    for my $gnm ( @gnms ) {
        my $fastaFile
            = -f "$fastaDir/$gnm.cds.gz"  ? "$fastaDir/$gnm.cds.gz"
            : -f "$fastaDir/$gnm.cds.bz2" ? "$fastaDir/$gnm.cds.bz2"
            : -f "$fastaDir/$gnm.cds"     ? "$fastaDir/$gnm.cds"
            : -f "$fastaDir/$gnm.ffn.gz"  ? "$fastaDir/$gnm.ffn.gz"
            : -f "$fastaDir/$gnm.ffn.bz2" ? "$fastaDir/$gnm.ffn.bz2"
            : -f "$fastaDir/$gnm.ffn"     ? "$fastaDir/$gnm.ffn"
            : "none";
        if( $fastaFile eq "none" ) {
            $cdsMistakes++;
            print "did not find CDS or FFN fasta file for $gnm\n";
        }
        else {
            my $openFile
                = $fastaFile =~ m{\.bz2} ? "bzip2 -qdc"
                : $fastaFile =~ m{\.gz}  ? "gzip -qdc"
                : "cat";
            my $id = "";
            open( my $FASTA,"-|","$openFile $fastaFile");
            while(<$FASTA>) {
                chomp;
                if( /^>(\S+)/ ) {
                    my $fullID = $1;
                    $id = $_ =~ m{protein_id=(\S+\.\d)\]} ? $1 : $fullID;
                    $id =~ s{lcl\|}{};
                    $id =~ s{_cds_\d+}{};
                }
                else {
                    $seq{"$id"} .= $_;
                }
            }
            close($FASTA);
            print "read fasta file:\n   $fastaFile\n";
        }
    }
    if( $cdsMistakes > 0 ) {
        return();
    }
    else {
        return(\%seq);
    }
}

sub getCodons {
    my $seq = $_[0];
    my @nucl = split(//,$seq);
    my @codons = ();
    my $codon = "";
    for my $nucl ( @nucl ) {
        $codon .= $nucl;
        if( length($codon) == 3 ) {
            push(@codons,$codon);
            $codon = "";
        }
    }
    return(@codons);
}

sub readPW {
    my $pwFile = $_[0];
    print "reading blast results:\n   $pwFile\n";
    my @lines2work = ();
    my %dNdSed     = ();
    open( my $ALN,"-|","bzip2 -dc $pwFile" );
  ALNLINE:
    while(<$ALN>) {
        next ALNLINE if( m{^#} );
        my(
            $query,$target,
            $evalue,$bitscore,$pident,
            $qstart,$qend,$qcov,
            $sstart,$send,$scov,
            $qalnseq,$salnseq
        ) = split;
        $query  = cleanID("$query");
        $target = cleanID("$target");
        next ALNLINE if( $query eq $target );
        next ALNLINE if( $dNdSed{"$query.$target"} > 0 );
        next ALNLINE if( $dNdSed{"$target.$query"} > 0 );
        $dNdSed{"$query.$target"}++;
        $dNdSed{"$target.$query"}++;
        my $line = join( "\t",
                         $query,$target,
                         $qstart,$qend,
                         $sstart,$send,
                         $qalnseq,$salnseq );
        push(@lines2work,$line);
    }
    close($ALN);
    my $count = @lines2work;
    if( $count > 0 ) {
        return($count,\@lines2work);
    }
    else {
        return($count);
    }
}

sub signalHandler {
    if( length($workDir) > 1 && -d "$workDir" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $workDir";
    }
    else {
        print "\n\tquitting $ownName\n";
    }
    die  "\tdone!\n\n";
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

sub inflateFile {
    my $origfile = $_[0];
    my( $opener ) = how2open($origfile);
    my $rootName = nakedName("$origfile");
    my $inflated = "$workDir/$rootName.faa";
    unless( -f "$inflated" ) {
        system("$opener $origfile > $inflated");
    }
    return("$inflated");
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
    my $inflated = "$workDir/$rootName.faa";
    unless( -f "$inflated" ) {
        system("$opener $origfile > $inflated");
    }
    return("$inflated");
}

sub progressLine {
    my($done,$toGo,$decim) = @_;
    my $columns = qx(tput cols);
    chomp($columns);
    if( $columns > 50
            && $toGo >= 10
            && -t STDOUT
            && -t STDIN ) {
        if( $decim > 2 ) {
            ### better to count than show percent
            my $buffer     = " " x ( length($toGo) - length($done) + 1 );
            my $counter    = "[$buffer" . $done . " ]";
            my $countSpace = length($counter);
            my $pbwidth    = $columns - ( $countSpace + 3 );
            my $nhashes    = int($done / $toGo * $pbwidth);
            printf("\r% -${pbwidth}s% ${countSpace}s",
                   '#' x $nhashes,$counter);
        }
        else {
            my $percent    = sprintf("%.${decim}f",(100 * $done / $toGo));
            my $maxPC      = sprintf("%.${decim}f",100);
            my $buffer     = " " x ( length($maxPC) - length($percent) + 1 );
            my $counter    = "[$buffer" . $percent . "% ]";
            my $countSpace = length($counter);
            my $pbwidth    = $columns - ( $countSpace + 3 );
            my $nhashes    = int($done / $toGo * $pbwidth);
            printf("\r% -${pbwidth}s% ${countSpace}s",
                   '#' x $nhashes,$counter);
        }
        if( $done == $toGo ) {
            print "\n";
        }
    }
}

sub cleanID {
    my $toClean = $_[0];
    $toClean =~ s{^\S+?\|(\S+)}{$1};
    $toClean =~ s{\|$}{};
    return("$toClean");
}

sub checkAlignSeqs {
    my $checkAlnFile = $_[0];
    my $itWorks = 0;
    open( my $ALNFILE,"-|","bzip2 -qdc $checkAlnFile" );
  ALNLINE:
    while(<$ALNFILE>) {
        if( m{^#} ) {
            s{^#\s*}{};
            my @testFields = split;
            my $countFields = @testFields;
            if( $testFields[-1] =~ m{^(tseq|sseq)$}i
                    && $testFields[-2] =~ m{^qseq$}i
                    && $countFields == 13 ) {
                $itWorks = 1;
            }
            last ALNLINE;
        }
    }
    close($ALNFILE);
    if( $itWorks > 0 ) {
        return($itWorks);
    }
    else {
        return();
    }
}

sub runPW {
    my($queryFile,$targetFile,$pwProg) = @_;
    print "running pairwise comparisons and RBH\n";
    my $command
        = qq(getRBH.pl -q $queryFile -t $targetFile )
        . qq( -p $pwProg -a T -x 2);
    my $pwOutput = qx($command);
}
