#!/usr/bin/perl
#zhijie.li@utoronto.ca
#2018 05 29
# takes two inputs: a fastafile name containing the query and a database name
#queries are in the fasta file
use strict;
use warnings;
use Cwd qw(abs_path);

#use YAML qw(Dump);
my $currentpath = abs_path '.';
#####read the input fasta, extract each record, send to species-sorting
open( my $input_fasta, "<", $ARGV[0] )
  or die "can't find input sequence file $ARGV[0]";
if ( $ARGV[0] !~ /.+\.fas(ta)*$/i ) {
    die "the first argument need .fas in filename ";
}

my $name;
my $name_cache;
my $desc;
my $sequence_ori;
my $count     = 0;
my $save_flag = 0;    #

my $singlename = '';
my %targets;
my $ser = sprintf( "%06d", $count );

while (<$input_fasta>) {
    chomp $_;
    if ( $_ =~ /^>(.+)$/ ) {
        $count++;
        $ser = sprintf( "%06d", $count );
        $targets{$ser}{'name'} = $1;
        $targets{$ser}{'seq'}  = '';

    }
    if ( $_ !~ /^>/ and $count > 0 ) {

        $targets{$ser}{'seq'} .= $_;
    }

}

###===end of fasta read cycle
close $input_fasta;

#open my $ofh, ">", "dump.yaml.txt";

#print $ofh Dump \%targets;
#close $ofh;
my $total = $count;

print "$total sequences read from file $ARGV[0] \n";

#prepare db
my $db = 'pdb_seqres.txt';
if ( defined $ARGV[1] ) {
    $db = $ARGV[1];
}

my $db_prefix = $db;
$db_prefix =~ s/\.(^[\.]+)$//;
my $db_unique = "$db_prefix.filtered.fas";
###blast treats da entries case insensitive, therefore must uniquify the entries
if ( !-d "blastdb" ) { mkdir "blastdb"; }
if ( !-e $db_unique ) { unique_fas( $db, $db_unique ); }

chdir "blastdb";

if ( !-e "..\\$db_unique" ) { die "can't find the file $db_unique"; }

if ( !-e "$db_prefix.psq" ) {
    system(
"makeblastdb -in \"..\\$db_unique\"  -dbtype prot -parse_seqids -out $db_prefix"
    );
}

#run blast
if ( -e "../$db_unique" and -e "$db_prefix.psq" ) {
    foreach my $k ( sort keys %targets ) {
        my $tar_fas = "$k.fas";

        open my $t_fh, ">", $tar_fas;
        print $t_fh ">$targets{$k}{'name'}\n$targets{$k}{'seq'}\n";
        close $t_fh;
        my $result = ` blastp -db $db_prefix -query $tar_fas `;
        test_hit( $result, "$targets{$k}{'name'}" );
        open my $logfh, '>>', "$k.blast_result.txt";
        print $logfh "$result\n";
        close $logfh;

    }
}
chdir $currentpath;

sub unique_fas {
    my ( $db, $db_unique ) = (@_);

    open my $da_fh, "<", $db;
    open my $ou_fh, ">", "$db_unique";
    my %ids;
    my $id = '';
    while (<$da_fh>) {
        my $line = $_;
        chomp $line;
        if ( $_ =~ /^>(.+)/ ) {
            if ( $_ =~ /^>([^\s\0]+)(.*)$/ ) {
                $id = lc $1;
                if ( !exists $ids{$id} ) {
                    $ids{$id} = 0;
                }
                else {
                    my $c     = 0;
                    my $newid = $id . '_' . $c;
                    while ( exists $ids{$newid} ) {
                        $c++;
                        $newid = $id . '_' . $c;
                        print "$newid\t";
                    }

                    #print "\n";
                    $ids{$newid} = 0;
                    $line = '>' . $newid . $2;
                }

            }
        }
        print $ou_fh "$line\n";

    }
    close $da_fh;
    close $ou_fh;
}

sub test_hit {
    my ( $po1result, $enzyme ) = @_;
    if ( $enzyme =~ /^(.+) OS/ ) { $enzyme = $1; }
    my @po1lines = split "\n", $po1result;

    my $evalue   = 1;
    my $bits     = 0;
    my $l        = '';
    my $db       = '';
    my $identity = 0;
    my $hit      = '';
    foreach my $i ( 0 .. ($#po1lines) ) {

        if ( ( index $po1lines[$i], '(Bits)  Value' ) > 0 ) {

            $l = $po1lines[ $i + 2 ];

#print "$po1lines[$i]\n[$l]\n";
#ref|XP_011062778.1|  PREDICTED: GDP-fucose protein O-fucosyltrans...  306     7e-102
            if ( $l =~ /^(.+)\s+([\d\.]+)\s+([e\d\-\.]+)\s*$/ ) {
                $hit    = $1;
                $bits   = $2;
                $evalue = $3;
            }
        }
        if (    $po1lines[$i] =~ /^>/
            and $po1lines[ $i + 4 ] =~ /Identities = \d+\/\d+ \((\d+)%\)/ )
        {
            $identity = $1;
            last;
        }

    }
    my $num_evalue = $evalue;
    if ( $evalue =~ /(\d)e\-([\d]+)/ ) { $num_evalue = $1 * 10**( -$2 ); }
    print ">$enzyme\nE-value:[$num_evalue]\t $hit\n";

    open my $fh1, ">>", "PDB_notfound.txt";
    open my $fh2, ">>", "PDB_high_homolog.txt";
    open my $fh3, ">>", "PDB_found.txt";
    open my $fh4, ">>", "Blast_result.txt";

    if ( $num_evalue > 0.01 ) {

# print "  $enzyme not found !!!!!!!!!         bits $bits   E-value $evalue ($num_evalue)\n";
        print $fh1 "identity: $identity%\t$bits\t$evalue\t$enzyme\t[$hit]\n";
    }
    elsif ( $num_evalue > ( 10**(-50) ) ) {

#print "  $enzyme domian containing??         bits $bits   E-value $evalue ($num_evalue)\n";
        print $fh2 "identity: $identity%\t$bits\t$evalue\t$enzyme\t[$hit]\n";
    }
    else {

#print "  $enzyme found                       bits $bits   E-value $evalue ($num_evalue)\n";
        print $fh3 "identity: $identity%\t$bits\t$evalue\t$enzyme\t[$hit]\n";
    }
    print $fh4 "identity: $identity%\t$bits\t$evalue\t$enzyme\t[$hit]\n";
    close $fh1;
    close $fh2;
    close $fh3;
    close $fh4;

}
