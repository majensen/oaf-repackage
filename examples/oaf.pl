#!/usr/bin/env perl
#
# Ornitine decarboxylase Antizyme Finder (OAF)
# Copyright 2006-2008 Bekaert M <mbekaert@gmail.com>
#
# This work is licensed under the Creative Commons Attribution-
# Noncommercial-Share Alike 3.0 License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use Bio::Seq;

use Bio::Tools::OAF;

#----------------------------------------
my $version = '1.0.6';
my $hmm     = getcwd() . '/oaz.hmm';
my $fs      = getcwd() . '/oaz.fs.hmm';
my $ref     = getcwd() . '/oaz.fasta';

#----------------------------------------------------------

sub hmm_disc {
    my ( $seq, $name, $id, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose ) = @_;

    my $chrom;
    if ( defined $name ) {
        $chrom = $name . '_' . $id;
    } elsif ( defined $seq->display_id ) {
        $chrom = $seq->display_id;
    } else {
        $chrom = 'noname_' . $id;
    }
    $seq->display_id($chrom);
    my $OAF_obj = Bio::Tools::OAF->new( $seq, $translate, $aug, $hmm, $fs );
    $OAF_obj->{'_verbose'} = 1 if ($verbose);
    $OAF_obj->{'_filter'} = $filter;
    if ( $OAF_obj->find($evalue) ) {
        $OAF_obj->show($format);
    } elsif ( $verbose && ( defined $OAF_obj->{'_hmmoaz'}[0] ) && ( $OAF_obj->{'_hmmoaz'}[0] < $OAF_obj->{'_evalue'} ) ) {
        $OAF_obj->show($format);
    }
}

sub fasta_filter {
    my ( $seq, $ref, $fasta, $organism, $id, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose ) = @_;

    my @hits = Bio::Tools::OAF->getHits( $seq, 1, $ref, $fasta );
    if ( $#hits >= 0 ) {
        for ( my $i = 0 ; $i <= $#hits ; $i++ ) {
            my ( $hitstart, $hitend, $hitstrand ) = split( /\|/, $hits[$i] );
            my $seqstart = $hitstart - 350;
            $seqstart = 1 if ( $seqstart < 1 );
            my $seqend = $hitend + 350;
            $seqend = $seq->length() if ( $seqend > $seq->length() );
            my $seq_rf = Bio::Seq->new(
                -seq      => $seq->subseq( $seqstart, $seqend ),
                -alphabet => 'dna',
                -id       => $seq->display_id . '|' . $id
            );
            my $OAF_obj = Bio::Tools::OAF->new( $seq_rf, $translate, $aug, $hmm, $fs );
            $OAF_obj->{'_verbose'} = 1 if ($verbose);
            $OAF_obj->{'_filter'} = $filter;

            if ( $OAF_obj->find( $evalue, $hitstrand, $seqstart, $seqend ) ) {
                $OAF_obj->show( $format, $organism );
            } elsif ( $verbose && ( defined $OAF_obj->{'_hmmoaz'}[0] ) && ( $OAF_obj->{'_hmmoaz'}[0] < $OAF_obj->{'_evalue'} ) ) {
                $OAF_obj->show( $format, $organism );
            }
        }
    }
}

#------------------------ Main ----------------------------
my ( $infile, $genbank, $translate, $name, $aug, $orthofile, $bank );
my ( $evalue, $verbose, $nofasta, $format, $filter, $remote, $out ) = ( 1e-40, 0, 0, 'fasta', 0, 0, 1 );

GetOptions(
    'sequence:s' => \$infile,
    'ortho:s'    => \$orthofile,
    'genbank:s'  => \$genbank,
    'a!'         => \$aug,
    'expect:f'   => \$evalue,
    'format:s'   => \$format,
    'name:s'     => \$name,
    'table:i'    => \$translate,
    'bank:s'     => \$bank,
    'd!'         => \$nofasta,
    'filter:i'   => \$filter,
    'r!'         => \$remote,
    'v!'         => \$verbose
);
$format = lc $format;

if (
    (
        ( ( defined $genbank ) && !( defined $orthofile ) && !( defined $infile ) ) 
        || (   !( defined $genbank )
            && !( defined $orthofile )
            && ( defined $infile )
            && ( -f $infile ) )
        || (   !( defined $genbank )
            && !( defined $infile )
            && ( defined $orthofile )
            && ( -f $orthofile ) )
    )
    && (
        !( defined $translate )
        || (   ( defined $translate )
            && ( $translate =~ /([1-9]|1[1-6]|2(1|2))/ ) )
    )
    && ( $filter >= 0 )
    && ( $filter <= 2 )
    && ( $format =~ /fasta|genbank|xml|raw/ )
  )
{
    $translate = 1     unless ( defined $translate );
    $bank      = 'est' unless ( defined $bank );
    print STDERR "\n..:: Ornitine decarboxylase Antizyme Finder (OAF) ::..\n> Standalone program version $version <\n\n"
      if ($verbose);
    if ( defined $infile ) {
        if ( my $inseq = Bio::SeqIO->new( '-file' => "<$infile", '-format' => 'fasta' ) ) {
            $out--;
            my $id = 0;
            Bio::Tools::OAF->show( 'xml-begin', $evalue, $translate )
              if ( $format eq 'xml' );
            while ( my $seq = $inseq->next_seq ) {
                $id++;
                if ( !$nofasta && ( length( $seq->seq() ) > 20000 ) ) {
                    fasta_filter( $seq, $ref, undef, undef, $id, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose );
                } else {
                    hmm_disc( $seq, $name, $id, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose );
                }
            }
            Bio::Tools::OAF->show('xml-end') if ( $format eq 'xml' );
        } else {
            print STDERR "FATAL: Incorrect file format.\n";
        }
    } elsif ( defined $genbank ) {
        my ( $translate, $organism, $seq ) = Bio::Tools::OAF->get_sequence($genbank);
        if ( defined $seq ) {
            $out--;
            Bio::Tools::OAF->show( 'xml-begin', $evalue, $translate )
              if ( $format eq 'xml' );
            if ( length( $seq->seq() ) > 20000 ) {
                fasta_filter( $seq, $ref, undef, $organism, 1, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose );
            } else {
                hmm_disc( $seq, $name, 1, $translate, $evalue, $format, $aug, $hmm, $fs, $filter, $verbose );
            }
            Bio::Tools::OAF->show('xml-end') if ( $format eq 'xml' );
        } else {
            $out += 2;
            print STDERR "FATAL: accession number not found!\n";
        }
    } elsif ( defined $orthofile ) {
        $out--;
        Bio::Tools::OAF->ortholog( $orthofile, $evalue, $format, $aug, $remote, $hmm, $fs, lc $bank, $filter, $verbose );
    } else {
        $out++;
    }
} else {
    print STDERR "\n..:: Ornitine decarboxylase Antizyme Finder (OAF) ::..\n> Standalone program version $version <\n\nFATAL: Incorrect arguments.\nUsage: oaf.pl [-options] mode\n\n Modes\n   --sequence=<sequence file>\n          Use a FASTA file (with one or more DNA/RNA sequences) as input.\n   --genbank=<accession>\n          Enable a search in GenBank. Provide an accession number rather\n          than a sequence file.\n   --ortho=<aa sequences file>\n          Start by searching for similar protein sequence(s) (FASTA file\n          format) in a database, and use those as queries for OAF.\n\n Options\n   -a\n         Force the use of alternative start codons, according to the\n         current genetic code table. Otherwise, ATG is the only allowed\n         initiation codon.\n   --expect\n         Set the E-value threshold. This setting specifies the statistical\n         significance threshold for reporting matches against database\n         sequences. The default value (1e-40).\n   --format\n         Available output format codes include 'fasta' (FASTA format);\n         'genbank' (Genbank format); 'raw' (raw sequence, no other\n         information); 'xml' (XML format). By default the FASTA format is\n         used.\n   --name\n         Overwrite the sequence name by the provided one. Otherwise the\n         program will use the sequence name from as input.\n   --table\n         Specify a genetic code table used for the translation of the query.\n         By default, for the 'genbank' and 'ortho' modes, the genetic code\n         table is used according to Genbank annotation. The 'sequence'\n         mode uses the Standard Code (1).\n\n Expert Options\n   --bank\n         Set the database used by the ortho[log] mode. For the remote\n         BLAST at NCBI (when -r is set) the available databases are\n         either 'refseq_rna' (NCBI Reference mRNA sequences database)\n         or 'est' (sequences from EST divisions). By default 'est' is\n         used.\n   -d\n         Turn off the pre-filter for long sequences and use user-provided\n         sequence directly for OAZ search. Normally, the sequences of\n         more than 20kb are pre-filtered by a local fasta search in order\n         to minimize the number of potential OAZ encoding sequences. This\n         option will usually increase the sensitivity of the search, if\n         the OAZ is VERY divergent, but that will dramatically increase\n         the time and memory used (so don't do it without good reason!)\n   --filter\n         Filter the sequences shown. Normally OAF.pl show all sequences\n         identified icluding putative/ambiguous which correspond to\n         ESTs containing truncated CDS, sequence errors or a missing +1\n         frameshift. By default there is no filter (0), you can filter\n         for the non-ambiguous hits only (1) or the ambiguous hits\n         only (2).\n   -r\n         Configure the program to use the remote BLAST server of the NCBI\n         rather than with a local BLAST.\n   -v\n         Print more possibly useful stuff, such as the individual scores\n         for each sequence.\n\n";
    $out++;
}

exit($out);
