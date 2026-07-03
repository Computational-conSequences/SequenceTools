#!/usr/bin/env perl

#########################################################
#						       	#
#	Author : Gabo Moreno-Hagelsieb         		#
#	       						#
#########################################################

use strict;
use Getopt::Long;
use Pod::Text;
use File::Temp qw( tempfile tempdir );
my $columns = qx(tput cols);
chomp($columns);
my $width = $columns >= 100 ? 80 : $columns - 3;
my $parser = Pod::Text->new (sentence => 0, width => $width, margin => 1);

my $ownName = $0;
$ownName =~ s{.*/}{};
my @famDBs = qw(
                   cog
                   cd
                   cdd
                   pfam
                   tigrfam
                   superfamily
                   VFDB
                   Toxins
           );
my $matchDom = join('|',@famDBs);

my @cddSets = qw(
                    cog
                    cd
                    cdd
            );
my $matchCDD = join("|",@cddSets);

my @progs = qw(
                  hmmscan
                  rpsblast
                  mmseqs
          );
my $mProgs = join("|",@progs);

####### defaults:
my @domFiles     = ();
my $famDB        = '';
my $outputFolder = 'cleanFams';
my $minCover     = 0.60;
my $maxOverlap   = 0.15;
my $appendAnn    = 'F';
my $reference    = '';

my $podUsage
    = qq(=pod\n\n)
    . qq(=head1 NAME\n\n)
    . qq($ownName - cleaning up results from scanFams.pl\n\n)
    . qq(=head1 SYNOPSIS\n\n)
    . qq($ownName -q [scanFamsFile] [options]\n\n)
    . qq(=head1 OPTIONS\n\n)
    . qq(=over\n\n)
    . qq(=item B<-q>\n\n)
    . qq(file(s) with family matches from scanFams.pl\n)
    . qq([e.g. GCF_000005845.cdd.mmseqs.bz2], required\n\n)
    . qq(=item B<-f>\n\n)
    . qq(family database, normally [$matchDom], make it explicit if not part\n)
    . qq(of the name of the file with match results\n)
    . qq([e.g. GCF_000005845.cdd.mmseqs.bz2]\n\n)
    . qq(=item B<-o>\n\n)
    . qq(output folder, default: $outputFolder\n\n)
    . qq(=item B<-c>\n\n)
    . qq(minimum coverage of family model, default 0.60\n\n)
    . qq(=item B<-v>\n\n)
    . qq(maximum overlap between matches, default 0.15\n\n)
    . qq(=item B<-r>\n\n)
    . qq(reference file with all scanned sequences [e.g. GCF_000005845.faa.gz]\n\n)
    . qq(=item B<-a>\n\n)
    . qq(append annotations (T|F). If 'T' the program will produce a file\n)
    . qq(with annotations appended (only works with $matchCDD),\n)
    . qq(default $appendAnn\n\n)
    . qq(=back\n\n)
    . qq(=head1 DESCRIPTION\n\n)
    . qq(B<This program> extracts results from hmmscan, rpsblast, or mmseqs\n)
    . qq(normally obtained using scanFams.pl\n\n)
    . qq(=cut\n\n)
    ;

GetOptions(
    "q=s{,}" => \@domFiles,
    "f=s"    => \$famDB,
    "o=s"    => \$outputFolder,
    "c=f"    => \$minCover,
    "v=f"    => \$maxOverlap,
    "r=s"    => \$reference,
    "a=s"    => \$appendAnn,
) or podhelp();
#pod2usage(1);

if ( !@domFiles ) {
    podhelp("I need files with resuls from scanFams.pl")
}

my @domFiles = do { my %seen; grep { !$seen{$_}++ } @domFiles };
my $countQueries = @domFiles;
my $foundQueries = 0;
my @missing      = ();
for my $domFile ( @domFiles ) {
    if( -f "$domFile" ) {
        $foundQueries++;
    }
    else {
        push(@missing,$domFile);
    }
}
if( $foundQueries == $countQueries ) {
    print "  found $countQueries scan files\n";
}
else {
    my $missing = @missing;
    my $error = "missing $missing matches file(s):\n".join("\n",@missing)."\n";
    die $error;
}

### the default is to not append annotations
my $appendAnn = $appendAnn =~ m{^(t|f)$}i ? uc($1) : "F";

if( length("$reference") > 1 ) {
    unless( -f "$reference" ) {
        podhelp("$reference file not found");
    }
}

#### genome DB:
my $genomeDB
    = exists $ENV{"GENOMEDB"} ? $ENV{"GENOMEDB"} : '.';

my $xfam
    = $famDB       =~ m{^($matchDom)$}i         ? uc($1)
    : $domFiles[0] =~ m{\.($matchDom)\.}i       ? uc($1)
    : $domFiles[0] =~ m{\.(\S+?)\.($mProgs)\.}i ? lc($1)
    : "none";
if( $xfam eq "none" ) {
    podhelp("I can only work with [$matchDom]");
}
print "  working with $xfam\n";
my $exten = lc($xfam);

#### COG specific
my $ncbiDir
    = -d "$genomeDB//DownLoad/ncbi" ? "$genomeDB//DownLoad/ncbi" : ".";
my $cddDir
    = -d "$ncbiDir/cdd" ? "$ncbiDir/cdd" : "none";
###### learn mini-functions for COGs and other CDDs
my %function  = ();
my %translate = ();
my %desc      = ();
my %priority  = ();
if( $xfam =~ m{^($matchCDD)$}i ) {
    my $cogsDir = $ncbiDir . "/COG/COG2024/";
    my $cogFunc = $cogsDir . "/data/cog-24.def.tab";
    print "  using COGs file:\n  $cogFunc\n";
    print "  learning COG functions\n";
    open( my $COGF,"<","$cogFunc" )
        or die "\tno COG functions file ($cogFunc)\n\n";
    while(<$COGF>) {
        next if( m{^#} );
        my($cog,$func,@descrip) = split;
        if( $cog =~ m{^COG} ) {
            $function{"$cog"} = $func;
        }
        else {
            print "\ttrouble with $_\n";
        }
    }
    close($COGF);
    ### learn how to translate families from CDD matches with rpsblast
    ###### learn CDD identifiers
    print "  learning $xfam identifiers from the CDD database\n";
    open( my $CDDF,"-|","gzip -qdc $cddDir/cddid_all.tbl.gz" )
        or die "\tno CDD functions file\n\n";
    while(<$CDDF>) {
        chomp;
        my($cdd,$collection,$secondary,$desc,$length) = split(/\t/,$_);
        $translate{"$cdd"}        = $collection;
        $translate{"$collection"} = $collection;
        if( $appendAnn eq "T" ) {
            $desc{"$collection"} = $desc;
        }
        if( $collection =~ m{^$xfam\d+}i ) {
            $priority{"$collection"}++;
        }
    }
    close($CDDF);
}

for my $domFile ( @domFiles ) {
    my $extractDir = $domFile =~ m{(\S+)/} ? $1 : "./";
    my $mainName
        = $domFile =~ m{$extractDir/(\S+)\.($matchDom)}i ? $1
        : $domFile =~ m{$extractDir/(\S+?)\.} ? $1
        : $domFile =~ m{$extractDir/(\S+)}   ? $1
        : $domFile =~ m{(\S+)} ? $1
        : "none";
    my $prog = $mainName =~ s{\.($mProgs)}{} ? $1 : "none";
    my $infile = $domFile;
    mkdir("$outputFolder") unless( -d "$outputFolder" );
    ###### open matches file
    ###### If working with COGs substitute identifiers for COGs with
    ###### functions
    if( $appendAnn eq "T" and $xfam =~ m{^($matchCDD)$}i ) {
        open( TRANSLATE,">","$outputFolder/$mainName.$exten.annot" );
    }
    my $translated = 0;
    my $original_format = 0;
    my %count
        = ( -f "$reference" ) ? learnReference("$reference")
        : ();
    my %checked = ();
    print "  working with $infile\n";
    my $format = 0;
    my ($cat,$trueFile) = figureCompression("$infile");
    open( my $DOMF,"-|","$cat $trueFile" );
  DOMLINES:
    while(<$DOMF>) {
        if( m{^#} ) {
            if( m{\t(Dcov|Fcov)} ) {
                $format = "clean";
                print "  the file's format is $format\n";
                next DOMLINES;
            }
            else {
                $format = $format eq "clean" ? $format : 3;
                print "  the file's format is $format (might be hmmscan)\n";
                next DOMLINES;
            }
        }
        my( $query,$dom_id,$score,
            $qStart,$qEnd,$qx,
            $domStart,$domEnd,$domcov
        ) = findFields($format,$_,$xfam,$prog);
        next DOMLINES if( $domcov < $minCover );
        my @toremember = ($qStart,$qEnd,$score,$domcov);
        if( $xfam eq "COG" ) {
            my $cog
                = length($translate{"$dom_id"}) > 1 ? $translate{"$dom_id"}
                : "NA";
            if( length($function{"$cog"}) > 0 ) {
                my $remember = join("\t",$cog,@toremember);
                push( @{$checked{"$query"}},$remember );
                $count{"$query"}++;
                $translated++;
                if( $appendAnn eq "T" ) {
                    my $newLine = $_;
                    my $append
                        = "\t" . $translate{"$dom_id"} . $function{"$cog"}
                        . "\t" . $desc{"$cog"};
                    $newLine =~ s{\n}{$append};
                    print TRANSLATE $newLine,"\n";
                }
            }
            else {
                # commented out because we know that some cogs
                # are missing in the new version for good reasons
                #print "\tproblems with $cog:\n$_";
            }
        }
        elsif( $xfam =~ m{^(CD|CDD)$} ) {
            if( length($translate{"$dom_id"}) > 2 ) {
                my $remember = join("\t",$translate{"$dom_id"},@toremember);
                push( @{$checked{"$query"}},$remember );
                $count{"$query"}++;
                $translated++;
                if( $appendAnn eq "T" ) {
                    my $newLine = $_;
                    my $id = $translate{"$dom_id"};
                    my $append
                        = "\t" . $id
                        . "\t" . $desc{"$id"};
                    $newLine =~ s{\n}{$append};
                    print TRANSLATE $newLine,"\n";
                }
            }
        }
        else {
            my $remember = join("\t",$dom_id,@toremember);
            push( @{ $checked{"$query"} },$remember );
            $count{"$query"}++;
            $translated++;
        }
    }
    close($DOMF);
    if( $appendAnn eq "T" ) {
        close(TRANSLATE);
    }
    if( $translated > 0 ) {
        if( $appendAnn eq "T" ) {
            system("bzip2 -9 $outputFolder/$mainName.$exten.annot");
        }
        open( my $TMPF,">","$outputFolder/$mainName.$exten" );
        print {$TMPF} join("\t","#Accs","${xfam}s(qstart:qend:domCover)"),"\n";
      CLEANQUERY:
        for my $query ( sort keys %count ) {
            if( $count{"$query"} == 0 ) {
                print {$TMPF} join("\t",$query,"NA"),"\n";
                next CLEANQUERY;
            }
            my @lines = @{ $checked{"$query"} };
            my %score = ();
            my %start = ();
            my %end   = ();
            for my $line ( @lines ) {
                my ($family,$qStart,$qEnd,$score,$domcov) = split(/\s+/,$line);
                #my $coords = join(":",$qStart,$qEnd);
                my $coords = join(":",$qStart,$qEnd,$domcov);
                my $fullMatch
                    = ( length($function{"$family"}) > 0 )
                    ? $family . $function{"$family"} . "($coords)"
                    : $family . "($coords)";
                $score{"$fullMatch"} = $score;
                $start{"$fullMatch"} = $qStart;
                $end{"$fullMatch"}   = $qEnd;
                $priority{"$fullMatch"} = $priority{"$family"};
            }
            my @acceptedCoords  = ();
            my @acceptedMatches = ();
          TESTMATCHES:
            for my $family (
                sort {
                    $priority{"$b"} <=> $priority{"$a"}
                        || $score{"$b"} <=> $score{"$a"}
                        || $start{"$a"} <=> $start{"$b"}
                        || $end{"$a"}   <=> $end{"$b"}
                        || $a <=> $b
                        || $a cmp $b
                    }
                keys %score ) {
                my $testStart = $start{"$family"};
                my $testEnd   = $end{"$family"};
                ##### allowing a bit of overlap:
                my $test_max_o
                    = sprintf("%.0f",
                              $maxOverlap * ($testEnd - $testStart + 1));
                for my $coords ( @acceptedCoords ) {
                    my($clearedStart,$clearedEnd) = split(/:/,$coords);
                    ##### allowing a bit of overlap:
                    my $cleared_max_o
                        = sprintf("%.0f",
                                  $maxOverlap *
                                  ($clearedEnd - $clearedStart + 1) );
                    my $acceptOverlap
                        = $cleared_max_o > $test_max_o ? $test_max_o
                        : $cleared_max_o;
                    #### start of overlap will be maximum start
                    my $coverStart
                        = $testStart >= $clearedStart ? $testStart
                        : $clearedStart;
                    #### end of overlap will be minimal end
                    my $coverEnd
                        = $testEnd <= $clearedEnd ? $testEnd : $clearedEnd;
                    my $overlap
                        = $testEnd < $clearedStart ? 0
                        : $clearedEnd < $testStart ? 0
                        : $coverEnd - $coverStart + 1;
                    next TESTMATCHES if( $overlap > $acceptOverlap );
                }
                ### add coords to accepted coords
                my $addCoords = join(":",$testStart,$testEnd);
                push(@acceptedCoords,$addCoords);
                ### add good non-overlapping matches to list
                push(@acceptedMatches,$family);
            }
            my @sortedMatches
                = sort { $start{"$a"} <=> $start{"$b"} } @acceptedMatches;
            my $printMatches = join(";",@sortedMatches);
            print {$TMPF} join("\t",$query,$printMatches),"\n";
        }
        close($TMPF);
        system("bzip2 -f -9 $outputFolder/$mainName.$exten");
    }
    else {
        print "   $domFile has no matching family results\n";
        unlink("$outputFolder/$mainName.$exten.annot");
    }
}

print "\tdone with $0\n\n";

sub findFields {
    my ($format,$line,$xfam,$prog) = @_;
    if( $format eq "clean" ) {
        my( $digest,$dom_id,$eval,$bitscore,
            $qStart,$qEnd,$qcov,
            $domStart,$domEnd,$domcov
        ) = split(/\s+/,$line);
        $digest =~ s{(lcl|ref)\|}{};
        $digest =~ s{\|$}{};
        $dom_id =~ s{gnl\|CDD\|}{};
        return( $digest,$dom_id,$bitscore,
                $qStart,$qEnd,$qcov,
                $domStart,$domEnd,$domcov);
    }
    elsif( $xfam =~ m{^($matchCDD)$}i
            || $prog =~ m{^(rpsblast|mmseqs)$} ) {
        my( $digest,$dom_id,$eval,$bitscore,
            $qStart,$qEnd,$q_ln,
            $domStart,$domEnd,$dom_len
        ) = split(/\s+/,$line);
        $digest =~ s{(lcl|ref)\|}{};
        $digest =~ s{\|$}{};
        $dom_id =~ s{gnl\|CDD\|}{};
        my $domcov = calcCoverage($domStart,$domEnd,$dom_len);
        return( $digest,$dom_id,$bitscore,
                $qStart,$qEnd,$q_ln,
                $domStart,$domEnd,$domcov);
    }
    else { ### it's an hmmscan file
        my( $dom_name,$dom_id,$dom_len,$digest,$caca1,$q_ln,
            $e1,$tscore,$tbias,$n1,$n2,
            $c_eval,$i_eval,$score,$bias,
            $domStart,$domEnd,
            $alnStart,$alnEnd,
            $qStart,$qEnd,$acc
        ) = split(/\s+/,$line);
        $digest =~ s{(lcl|ref)\|}{};
        $digest =~ s{\|$}{};
        my $dom
            = $dom_id eq "-" ? $dom_name
            : $dom_id;
        my $domcov = calcCoverage($domStart,$domEnd,$dom_len);
        return( $digest,$dom,$score,
                $qStart,$qEnd,$q_ln,
                $domStart,$domEnd,$domcov);
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

sub learnReference {
    my $inref = $_[0];
    my %ids   = ();
    if( -f "$inref" ) {
        my $openref
            = $inref =~ m{\.bz2}    ? "bzip2 -qdc $inref"
            : $inref =~ m{\.(gz|Z)} ? "gzip -qdc $inref"
            : "cat $inref";
        open( my $INREF,"-|","$openref" );
        while(<$INREF>) {
            if( m{^>(\S+)} ) {
                $ids{"$1"} = 0;
            }
        }
        close($INREF);
        my $countids = keys %ids;
        if( $countids > 0 ) {
            return(%ids);
        }
        else {
            die "reference file $inref had no sequence identifiers\n\n";
        }
    }
    else {
        die "I did not find the reference sequence file: $inref\n\n";
    }
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
        = $fullName =~ m{\.(gz|Z)$} ? "gzip  -qdc"
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
