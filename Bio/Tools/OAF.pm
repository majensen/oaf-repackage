## $Id: OAF.pm,v 1.0.6 2008/01/25 02:53:11 Europe/Dublin $
#
# Ornitine decarboxylase Antizyme Finder (OAF)
# Copyright 2006-2008 Bekaert M <mbekaert@gmail.com>
#
# This work is licensed under the Creative Commons Attribution-
# Noncommercial-Share Alike 3.0 License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::OAF - Ornithine decarboxylase Antizyme Finder (bioperl module)

=head1 SYNOPSIS

Take a sequence object from, say, an inputstream, and find an eukaryotic
ornithine decarboxylase antizyme. HMM profiles are used in order to
identify location, frame and orientation of such gene.

Creating the OAF object, eg:

  my $inputstream = Bio::SeqIO->new( -file => 'seqfile', -format => 'fasta' );
  my $seqobj = $inputstream->next_seq();
  my $OAF_obj = Bio::Tools::OAF->new( $seqobj );

Obtain an array holding the start point, the stop point, the DNA sequence
and amino-acid sequence, eg:

  my $array_ref = $OAF_obj->{'_result'} if ( $OAF_obj->find() );


Display result in genbank format, eg:

  $OAF_obj->show( 'genbank' );

=head1 DESCRIPTION

Bio::Tools::OAF is a featherweight object for identifying ornithine
decarboxylase antizymes. The sequences should not be longer than 20kb. The
returned array include location, sequence and statistic for the putative
ornithine decarboxylase antizyme gene. Fully functional with mRNA and EST
sequence, no intron allowed.

See Synopsis above for the object creation code.

=head1 DRIVER SCRIPT

  #!/usr/bin/perl -w
  use strict;

  use Bio::Seq;
  use Bio::Tools::OAF;

  my $inseq = Bio::SeqIO->new( '-file' => "< $yourfile", -format => 'fasta' );
  while (my $seq = $inseq->next_seq) {
    my $OAF_obj = Bio::Tools::OAF->new( $seq, 11 );
    if ( $OAF_obj->find() ) {
      $OAZ_obj->show( 'genbank' );
    } else {
      print "  no hit!\n";
    }
  }

=head1 REQUIREMENTS

To use this module you may need:
 * Bioperl (L<http://www.bioperl.org/>) modules,
 * HMMER distribution (L<http://hmmer.janelia.org/>),
 * Blast (L<http://www.ncbi.nlm.nih.gov/BLAST/download.shtml>) and
 * FASTA distribution (L<ftp://ftp.ebi.ac.uk/pub/software/unix/fasta/>).

=head1 FEEDBACK

User feedback is an integral part of the evolution of this modules. Send
your comments and suggestions preferably to author.

=head1 AUTHOR

B<Michael Bekaert> (mbekaert@gmail.com)

Address:
     School of Biology & Environmental Science
     University College Dublin
     Belfield, Dublin 4
     Dublin
     Ireland

=head1 SEE ALSO

perl(1), manual.html, bioperl web site

=head1 LICENSE

Copyright 2006-2008 - Michael Bekaert

This work is licensed under the Creative Commons Attribution-
Noncommercial-Share Alike 3.0 License. To view a copy of this
license, visit L<http://creativecommons.org/licenses/by-nc-sa/3.0/>
or send a letter to Creative Commons, 543 Howard Street, 5th
Floor, San Francisco, California, 94105, USA

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

package Bio::Tools::OAF;

use strict;
use vars qw(@ISA $VERSION);
use File::Temp qw/tempfile/;
use Bio::DB::GenBank;
use Bio::Tools::HMMER::Results;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use Bio::PrimarySeq;
use Bio::PrimarySeqI;
use Bio::Tools::CodonTable;

use Bio::Tools::BlastLocalTools;
use Bio::Tools::BlastTools;

$VERSION = '1.0.6';
@ISA     = qw(Bio::Root::Root Bio::Root::IO);

# Default path
my $PATH_REF = './oaz.fasta';
my $PATH_HMM = './oaz.hmm';
my $PATH_FS  = './oaz.fs.hmm';
my $PATH_TMP = '/tmp';

=head2 _findexec

 Title   : _findexec
 Usage   : my $path = $OAF_obj->_findexec( $exec );
 Function: Find an executable file in the $PATH.
 Returns : The full path to the executable.
 Args    : $exec (mandatory) executable to be find.

=cut                       

sub _findexec {
    my ( $self, @args ) = @_;
    my $exec = shift(@args);

    foreach my $p ( split( /:/, $ENV{'PATH'} ) ) {
        return "$p/$exec" if ( -x "$p/$exec" );
    }
}

=head2 new

 Title   : new
 Usage   : my $nb = Bio::Tools::OAF->new( $seqobj, $table, $aug, $hmm, $fs );
 Function: Initialize OAF object.
 Returns : An OAF object.
 Args    : $seqobj (mandatory) PrimarySeqI object (dna or rna),
           $table (optional) translation table/genetic code number,
              the default value is 1,
           $aug (optional) use other start codon than AUG (default 0),
           $hmm (optional) path to hmm profiles by default OAF looks at
             ./oaz.hmm
           $fs (optional) path to frameshift hmm profiles by default OAF
             looks at ./oaz.fs.hmm.

=cut

sub new {
    my ( $self, @args ) = @_;
    my ( $seqobj, $table, $aug, $hmm, $fs ) = @args;
    $self = {};

    $self->{'_aug'}  = ( ( defined $aug && int($aug) == 1 ) ? 0                : 1 );
    $self->{'_hmm'}  = ( ( defined $hmm )                   ? $hmm             : $PATH_HMM );
    $self->{'_fs'}   = ( ( defined $fs )                    ? $fs              : $PATH_FS );
    $self->{'_pfam'} = ( ( defined $ENV{'HMMERDIR'} )       ? $ENV{'HMMERDIR'} : Bio::Tools::OAF->_findexec('hmmpfam') );
    $seqobj->throw("die in _initialize, OAF.pm works only on PrimarySeqI objects\n")
      unless ( $seqobj->isa("Bio::PrimarySeqI") );
    $seqobj->throw("die in _initialize, OAF.pm works only on DNA sequences\n")
      if ( $seqobj->alphabet eq 'protein' );
    $seqobj->throw("die in _initialize, OAF.pm works only on DNA sequences < 20kb\n")
      if ( length( $seqobj->seq() ) > 20000 );
    $seqobj->throw("die in _initialize, hmmpfam not found\n")
      unless ( -x $self->{'_pfam'} );
    $seqobj->throw( 'die in _initialize, hmm profile not found at ' . $self->{'_hmm'} . "\n" )
      unless ( -f $self->{'_hmm'} );
    $seqobj->throw( 'die in _initialize, hmm/fs profile not found at ' . $self->{'_fs'} . "\n" )
      unless ( -f $self->{'_fs'} );
    my $chrom = ( ( defined $seqobj->display_id ) ? $seqobj->display_id : 'noname' );
    $seqobj = uc $seqobj->seq();
    $seqobj =~ tr/U/T/;
    $self->{'_seqref'} = Bio::Seq->new( -seq => $seqobj, -alphabet => 'dna', -id => $chrom );
    $self->{'_table'} = ( ( defined $table ) ? $table : 1 );
    $self->{'_filter'} = 0;
    bless($self);
    return $self;
}

=head2 find

 Title   : find
 Usage   : my $bool = $OAF_obj->find( $evalue, $strand, $start, $end );
 Function: Identify an ornithine decarboxylase antizyme protein.
 Returns : boolean.
 Args    : $evalue (optional) set the E-value (expected) threshold.
             Default is 1e-40,
           $strand(optional) strand where search should be done (1 direct,
             -1 reverse or 0 both). Default is 0,
           $start (optional) coordinate of the first nucleotide. Useful
             for coordinate calculation when first is not 1. Default is 1,
           $end (optional) coordinate of the last nucleotide. Default is
             the sequence length.

=cut

sub find {
    my ( $self, @args ) = @_;
    my ( $evalue, $strand, $start, $end ) = @args;

    $self->{'_evalue'} = ( ( defined $evalue ) && ( $evalue > 0 ) ) ? $evalue : 1e-40;
    $strand = 0
      unless ( ( defined $strand ) && ( $strand == 1 || $strand == -1 ) );
    $start = 1 unless ( ( defined $start ) && $start > 1 );
    $end = length( $self->{'_seqref'}->seq() )
      unless ( ( defined $end ) && $end > 1 && $end > $start );
    my $status = $self->_what_oaz($strand);
    if ( ($status) && ( $self->{'_hmmoaz'}[0] < $self->{'_evalue'} ) ) {
        print STDERR "\n> ", ( substr $self->{'_hmmoaz'}[1], 0, 4 ), ' found (', $self->{'_hmmoaz'}[0], ') in frame ', $self->{'_hmmoaz'}[2], ( ( ( defined $self->{'_hmmfs'}[0] ) && ( $self->{'_hmmoaz'}[2] != $self->{'_hmmfs'}[2] ) ) ? ( ' and first ORF in frame ' . $self->{'_hmmfs'}[2] ) : '' ), "\n"
          if ( $self->{'_verbose'} );
        return $self->_find_orf( $self->{'_hmmoaz'}[6], $start, $end );
    }
    return 0;
}

=head2 _what_oaz

 Title   : _what_oaz
 Usage   : my $bool = $OAF_obj->_what_oaz( $strand );
 Function: Use HMM profiles to identify an ornithine decarboxylase antizyme protein.
 Returns : boolean.
 Args    : $strand(optional) strand where search should be done (1 direct, -1 reverse or 0 both). Default is 0.

=cut

sub _what_oaz {
    my ( $self, @args ) = @_;
    my $strand = shift(@args);

    my $seq = $self->{'_seqref'};
    $seq = $seq->revcom if ( $strand < 0 );

    my ( $TMP, $filename ) = tempfile( DIR => $PATH_TMP, UNLINK => 1 );
    for ( my $i = 0 ; $i < 3 ; $i++ ) {
        print $TMP ">$i\n", $seq->translate( undef, undef, $i, $self->{'_table'} )->seq(), "\n";
    }
    if ( $strand == 0 ) {
        $seq = $seq->revcom;
        for ( my $i = 0 ; $i < 3 ; $i++ ) {
            print $TMP ">$i-\n", $seq->translate( undef, undef, $i, $self->{'_table'} )->seq(), "\n";
        }
    }
    close $TMP;
    system( $self->{'_pfam'} . ' ' . $self->{'_hmm'} . ' ' . $filename . ' > ' . $filename . '.report' );
    my $hmmer = new Bio::Tools::HMMER::Results(
        -file => $filename . '.report',
        -type => 'hmmpfam'
    );
    foreach my $domain ( $hmmer->each_Domain ) {
        $self->{'_hmmoaz'} = (
            [
                $domain->evalue,
                $domain->hmmname,
                substr( $domain->seq_id, 0, 1 ),
                $domain->bits,
                $domain->start,
                $domain->end,
                (
                    ( $strand != 0 )
                    ? $strand
                    : ( ( length( $domain->seq_id ) > 1 ) ? -1 : 1 )
                )
            ]
          )
          if ( ( ( !defined $self->{'_hmmoaz'}[0] ) || $domain->evalue < $self->{'_hmmoaz'}[0] )
            && ( ( substr $domain->hmmname, -5, 5 ) ne '_ORF0' ) );
    }
    foreach my $domain ( $hmmer->each_Domain ) {
        $self->{'_hmmfs'} = (
            [
                $domain->evalue,
                $domain->hmmname,
                substr( $domain->seq_id, 0, 1 ),
                $domain->bits,
                $domain->start,
                $domain->end,
                (
                    ( $strand != 0 )
                    ? $strand
                    : ( ( length( $domain->seq_id ) > 1 ) ? -1 : 1 )
                )
            ]
          )
          if (
            ( defined $self->{'_hmmoaz'}[0] )
            && ( ( !defined $self->{'_hmmfs'}[0] )
                || $domain->evalue < $self->{'_hmmfs'}[0] )
            && ( ( substr $domain->hmmname, -5, 5 ) eq '_ORF0' )
            && ( ( substr $domain->hmmname, 3, 1 ) eq ( substr $self->{'_hmmoaz'}[1], 3, 1 ) )
          );
    }
    unlink($filename);
    unlink( $filename . '.report' );
    return ( defined $self->{'_hmmoaz'}[0] ) ? 1 : 0;
}

=head2 _what_fs

 Title   : _what_fs
 Usage   : my $evalue = $OAF_obj->_what_fs( $site );
 Function: Use an HMM profile to qualify the -1 frameshift site.
 Returns : double (e-value).
 Args    : $site (mandatory) the frameshift site sequence.

=cut

sub _what_fs {
    my ( $self, @args ) = @_;
    my $site = shift(@args);
    my $evalue;

    my ( $TMP, $filename ) = tempfile( DIR => $PATH_TMP, UNLINK => 1 );
    print $TMP ">fs\n", $site, "\n";
    close $TMP;
    system( $self->{'_pfam'} . ' -n ' . $self->{'_fs'} . ' ' . $filename . ' > ' . $filename . '.report' );
    my $hmmer = new Bio::Tools::HMMER::Results(
        -file => $filename . '.report',
        -type => 'hmmpfam'
    );

    foreach my $domain ( $hmmer->each_Domain ) {
        $evalue = $domain->evalue if ( $domain->hmmname eq 'fs' );
    }
    unlink($filename);
    unlink( $filename . '.report' );
    return $evalue;
}

=head2 _find_orf

 Title   : _find_orf
 Usage   : my $bool = $OAF_obj->_find_orf( $strand, $start, $end );
 Function: Retrieve the ornithine decarboxylase antizyme ORF.
 Returns : boolean.
 Args    : $strand (mandatory) Strand where OAZ have been found
           (1 direct or -1 reverse),
           $start (mandatory) coordinate of the first nucleotide,
           $end (mandatory) coordinate of the last nucleotide.

=cut

sub _find_orf {
    my ( $self, @args ) = @_;
    my ( $strand, $start, $end ) = @args;

    my %name = ( 'L', 'AZL', 'S', 'AZS', 'A', 'AZ', 'C', 'AZ', '2', 'OAZ2', 'F', 'AZ', '1', 'OAZ1', '3', 'OAZ3', 'M', 'AZ', 'N', 'AZ' );
    my ( $position1, $position2 );
    my ( $begin, $stop ) = $self->_translation();
    my $seq = $self->{'_seqref'};
    $seq = $seq->revcom if ( $strand < 0 );
    $seq = $seq->seq();
    if (   ( defined $self->{'_hmmoaz'} )
        && ( defined $self->{'_hmmfs'}[0] )
        && ( $self->{'_hmmoaz'}[2] != $self->{'_hmmfs'}[2] ) )
    {

        #OAZ with FS
        my $fs = 1 + ( $self->{'_hmmoaz'}[2] == ( $self->{'_hmmfs'}[2] - 1 ) % 3 );
        $self->{'_minus'} = 1 if ( $fs == 2 );
        my $mydna = substr( $seq, 0, ( $self->{'_hmmfs'}[4] * 3 + abs( $self->{'_hmmfs'}[2] ) ) );
        if ( $mydna =~ m/(($begin)((?!($stop|$begin))(.{3}))*?)$/o ) {
            $position1 = length($`);
        } else {
            if ( $mydna =~ m/(($stop)((?!($stop))(.{3}))*)$/o ) {
                $position1 = length($`) + 3;
                $self->{'_noaug'} = 1;
            } elsif ( $mydna =~ m/^(.{0,2})((.{3})*)$/o ) {
                $self->{'_no5'} = 1 if ( $self->{'_hmmfs'}[4] < 10 );
                $position1 = length($1);
                $self->{'_noaug'} = 1;
            }
        }
        if ( defined $position1 ) {
            my $mydna = substr( $seq, $position1 );
            if (
                (
                    ( $fs == 1 ) && ( ( $mydna =~ m/^(.{3})(((?!($stop))(.{3}))*)($stop)(.{1}((?!($stop))(.{3})){75,})($stop)/o )
                        || ( $mydna =~ m/^(.{3})(((?!($stop))(.{3}))*)($stop)(.{1}((?!($stop))(.{3})){75,})(.{0,2})$/o ) )
                )
                || (
                    ( $fs == 2 )
                    && (   ( $mydna =~ m/^(.{3})(((?!($stop))(.{3}))*)($stop)(.{2}((?!($stop))(.{3})){75,})($stop)/o )
                        || ( $mydna =~ m/^(.{3})(((?!($stop))(.{3}))*)($stop)(.{2}((?!($stop))(.{3})){75,})(.{0,2})$/o ) )
                )
              )
            {
                my ( $position3, $coord );
                my $stopcodon = ( ( length($11) == 3 ) ? $11 : '' );
                $self->{'_no3'} = 1 if ( length($11) < 3 );
                if ( $strand > 0 ) {
                    $position2 = $start + $position1 + length( $1 . $2 . $6 . $7 . $stopcodon ) - 1;
                    $position3 = $start + $position1 + length( $1 . $2 );
                    $position1 = $start + $position1;
                    $coord     = 'join(' . ( ( defined $self->{'_no5'} ) ? '<' : '' ) . $position1 . '..' . ( $position3 - 1 ) . ',' . ( $position3 + 3 + ( -2 * $fs ) ) . '..' . $position2 . ( ( defined $self->{'_no3'} ) ? '>' : '' ) . ')';
                } else {
                    $position2 = $end - $position1;
                    $position3 = $end - $position1 - length( $1 . $2 );
                    $position1 = $end - $position1 - length( $1 . $2 . $6 . $7 . $stopcodon ) + 1;
                    $coord     = 'complement(join(' . ( ( defined $self->{'_no3'} ) ? '<' : '' ) . $position1 . '..' . ( $position3 - 1 ) . ',' . ( $position3 + 3 + ( -2 * $fs ) ) . '..' . $position2 . ( ( defined $self->{'_no5'} ) ? '>' : '' ) . '))';
                }
                $self->{'_result'} = (
                    [
                        $coord,
                        $position1,
                        $position2,
                        $strand,
                        'Ornithine decarboxylase antizyme ' . $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) },
                        $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) },
                        ( ( $fs == 2 ) ? 'Putative -1 frameshift; ' : '' )
                          . (
                            ( defined $self->{'_no5'} ) ? 'The sequence seems incomplete, 5\' of the CDS is missing; '
                            : ''
                          )
                          . (
                            ( defined $self->{'_no3'} ) ? 'The sequence seems incomplete, 3\' of the CDS is missing; '
                            : ''
                          )
                          . (
                            (
                                ( defined $self->{'_noaug'} ) ? 'The position of the initiation codon is not identified; '
                                : ''
                            )
                            . $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) }
                              . ' ORF0: '
                              . $self->{'_hmmfs'}[0]
                              . ' ORF1: '
                              . $self->{'_hmmoaz'}[0]
                          ),
                        $1 . $2 . $6 . $7 . $stopcodon,
                        Bio::Seq->new(
                            -seq => $1 . $2 . ( ( $fs == 1 ) ? substr( $6, 1 ) : substr( $2, -1 ) . $6 ) . $7,
                            alphabet => 'dna'
                          )->translate( undef, undef, undef, $self->{'_table'} )->seq()
                    ]
                );
            }
        }
    } elsif ( defined $self->{'_hmmoaz'} ) {

        #OAZ without FS
        $self->{'_no5'} = 1 if ( $self->{'_hmmoaz'}[4] < 10 );
        my $mydna = substr( $seq, 0, ( $self->{'_hmmoaz'}[4] * 3 + abs( $self->{'_hmmoaz'}[2] ) + 9 ) );
        if ( $mydna =~ m/(($begin)((?!($stop|$begin))(.{3}))*?)$/o ) {
            $position1 = length($`);
        } else {
            if ( $mydna =~ m/(($stop)((?!($stop))(.{3}))*)$/o ) {
                $position1 = length($`) + 3;
                $self->{'_noaug'} = 1;
            } elsif ( $mydna =~ m/^(.{0,2})((.{3})*)$/o ) {
                $position1 = length($1) + 1;
                $self->{'_noaug'} = 1;
            }
        }
        if ( defined $position1 ) {
            my $mydna = substr( $seq, $position1 );
            my ( $coord, $dna, $stopcodon );
            if ( $mydna =~ m/^(.{3})(((?!($stop))(.{3})){75,})($stop)/o ) {
                $dna       = $1 . $2;
                $stopcodon = $6;
            } elsif ( $seq =~ m/($stop)(((?!($stop))(.{3})){250,})($stop)/o ) {
                $self->{'_noaug'} = 1;
                $position1        = length($`) + 3;
                $dna              = $2;
                $stopcodon        = $6;
            } elsif ( $mydna =~ m/^(.{3})(((?!($stop))(.{3})){75,})(.{0,2})$/o ) {
                $self->{'_no3'} = 1;
                $dna            = $2;
                $stopcodon      = '';
            }
            if ( defined $dna ) {
                if ( $strand > 0 ) {
                    $position2 = $start + $position1 + length( $dna . $stopcodon ) - 1;
                    $position1 = $start + $position1;
                    $coord     = ( ( defined $self->{'_no5'} ) ? '<' : '' ) . $position1 . '..' . $position2 . ( ( defined $self->{'_no3'} ) ? '>' : '' );
                } else {
                    $position2 = $end - $position1;
                    $position1 = $end - $position1 - length( $dna . $stopcodon ) + 1;
                    $coord     = 'complement(' . ( ( defined $self->{'_no3'} ) ? '<' : '' ) . $position1 . '..' . $position2 . ( ( defined $self->{'_no5'} ) ? '>' : '' ) . ')';
                }
                $self->{'_result'} = (
                    [
                        $coord,
                        $position1,
                        $position2,
                        $strand,
                        'Ornithine decarboxylase antizyme ' . $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) },
                        $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) },
                        'No frameshift identified; '
                          . (
                            ( defined $self->{'_no5'} ) ? 'The sequence seems incomplete, 5\' of the CDS is missing; '
                            : ''
                          )
                          . (
                            ( defined $self->{'_no3'} ) ? 'The sequence seems incomplete, 3\' of the CDS is missing; '
                            : ''
                          )
                          . ( ( !( defined $self->{'_no5'} ) && !( defined $self->{'_no3'} ) ) ? 'The sequence seems erroneous; ' : '' )
                          . (
                            (
                                ( defined $self->{'_noaug'} ) ? 'The position of the initiation codon is not identified; '
                                : ''
                            )
                            . $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) }
                              . ' ORF: '
                              . $self->{'_hmmoaz'}[0]
                          ),
                        $dna . $stopcodon,
                        Bio::Seq->new( -seq => $dna, alphabet => 'dna' )->translate( undef, undef, undef, $self->{'_table'} )->seq()
                    ]
                ) if ( defined $coord );
            }
        }
    }
    if ( ( defined $self->{'_hmmoaz'} ) && !( defined $self->{'_result'} ) ) {
        $self->{'_result'} = ( [ ( $strand > 0 ) ? $start . '..' . $end : 'complement(' . $start . '..' . $end . ')', $start, $end, $strand, '', $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) }, 'The sequence seems incomplete, the CDS is missing; ' . 'Putative ' . $name{ ( substr $self->{'_hmmoaz'}[1], 3, 1 ) } . ' ORF' . ( ( ( defined $self->{'_hmmfs'}[0] ) && ( $self->{'_hmmoaz'}[2] != $self->{'_hmmfs'}[2] ) ) ? '0: ' . $self->{'_hmmfs'}[0] . ' ORF1: ' . $self->{'_hmmoaz'}[0] : ': ' . $self->{'_hmmoaz'}[0] ), $seq ] );

    }
    return ( defined $position2 ) ? 1 : 0;
}

=head2 getHits

 Title   : getHits
 Usage   : my @hits = Bio::Tools::OAF->getHits( $evalue, $seq, $ref );
 Function: Quick localization of OAZs (use FASTA).
 Returns : Array of hits start/stop positions.
 Args    : $evalue (mandatory) det the E-value threshold,
           $seq (mandatory) primarySeqI object (dna or rna),
           $ref (optional) path to fasta reference file, by default OAF
             look at ./OAZ.fasta.

=cut

sub getHits {
    my ( $self, @args ) = @_;
    my ( $seq, $evalue, $ref ) = @args;
    $ref = ( ( defined $ref ) ? $ref : $PATH_REF );
    $evalue = ( ( defined $evalue ) && ( $evalue > 0 ) ) ? $evalue : 1;

    my ( $TMP, $filename ) = tempfile( DIR => $PATH_TMP, UNLINK => 1 );
    my $fasta = ( ( defined $ENV{'FASTADIR'} ) ? $ENV{'FASTADIR'} : $self->_findexec('tfasta34') );
    my @hits;
    if ( ( defined $seq ) && ( -x $fasta ) ) {
        print $TMP ">query\n", $seq->seq(), "\n\n";
        close($TMP);
        system( $fasta . ' -Q -b 1 -d 1 -H ' . $ref . ' ' . $filename . '> ' . $filename . '.report' );
        my $in = new Bio::SearchIO(
            -format => 'fasta',
            -file   => $filename . '.report'
        );
        eval {
            while ( my $result = $in->next_result )
            {
                while ( my $hit = $result->next_hit ) {
                    while ( my $hsp = $hit->next_hsp ) {
                        push( @hits, ( $hsp->start('hit') . '|' . $hsp->end('hit') . '|' . $hit->strand('query') ) ) if ( $hsp->evalue <= $evalue );
                    }
                }
            }
        };
        unlink( $filename . '.report' );
        unlink($filename);
    }
    if ( $#hits >= 0 ) {
        @hits = sort @hits;
        my @hits2 = ();
        for ( my $i = 0 ; $i < $#hits ; $i++ ) {
            my ( $hitstart, $hitend, $hitstrand ) = split( /\|/, $hits[$i] );
            my $loop = 0;
            do {
                $loop = 0;
                my ( $hitstart_n, $hitend_n, $hitstrand_n ) =
                  split( /\|/, $hits[ $i + 1 ] )
                  if ( defined $hits[ $i + 1 ] );

                #frameshift or echo ?
                if (
                    ( defined $hits[ $i + 1 ] )
                    && (   ( abs( $hitstart_n - $hitstart ) < 1500 )
                        || ( abs( $hitend_n - $hitend ) < 1500 )
                        || abs( $hitstart_n - $hitend ) < 1500 )
                    && ( $hitstrand == $hitstrand_n )
                  )
                {
                    $i++;
                    $hitend   = $hitend_n   if ( $hitend < $hitend_n );
                    $hitstart = $hitstart_n if ( $hitstart > $hitstart_n );
                    $loop     = 1;
                }
            } until ( $loop == 0 );
            push( @hits2, ( $hitstart . '|' . $hitend . '|' . $hitstrand ) );
        }
        @hits = @hits2;
    }
    return @hits;
}

=head2 show

 Title   : show
 Usage   : $OAF_obj->show( $outstyle );
 Function: Print result in various style.
 Returns : none.
 Args    : $outstyle (mandatory) 'fasta', 'genbank', 'raw' or 'xml' style.

=cut

sub show {
    my ( $self, @args ) = @_;
    my $out = shift(@args);

    if ( $out eq 'xml-begin' ) {
        print "<oafxml version=\"1.2\">\n";
        print " <analysis>\n";
        print "  <program>\n";
        print "   <prog-name>OAF.pm</prog-name>\n";
        print "   <prog-version>$VERSION</prog-version>\n";
        print "  </program>\n";
        print "  <date>\n";
        print '   <day>', (gmtime)[3], "</day>\n";
        print '   <month>', (gmtime)[4] + 1,    "</month>\n";
        print '   <year>',  (gmtime)[5] + 1900, "</year>\n";
        print "  </date>\n";
        print "  <parameter>\n";
        print '   <evalue>', shift(@args), "</evalue>\n";
        print '   <table>',  shift(@args), "</table>\n";
        print "  </parameter>\n";
        print " </analysis>\n";
    } elsif ( $out eq 'xml-end' ) {
        print "</oafxml>\n";
    } elsif ( ( defined $self->{'_result'} ) && ( defined $self->{'_result'}[8] ) && ( $self->{'_filter'} == 0 || ( ( $self->{'_filter'} == 1 ) && ( $self->{'_result'}[6] =~ /ORF1/ ) && !( defined $self->{'_minus'} ) && !( defined $self->{'_no3'} ) && !( defined $self->{'_no5'} ) ) || ( ( $self->{'_filter'} == 2 ) && ( !( $self->{'_result'}[6] =~ /ORF1/ ) || ( defined $self->{'_minus'} ) || ( defined $self->{'_no3'} ) || ( defined $self->{'_no5'} ) ) ) ) ) {
        if ( $out eq 'genbank' ) {
            print '     gene            ',
              (
                  ( $self->{'_result'}[3] < 0 )
                ? ( 'complement(' . ( ( defined $self->{'_no3'} ) ? '<' : '' ) . $self->{'_result'}[1] . '..' . $self->{'_result'}[2] . ( ( defined $self->{'_no5'} ) ? '>' : '' ) . ')' )
                : ( ( ( defined $self->{'_no5'} ) ? '<' : '' ) . $self->{'_result'}[1] . '..' . $self->{'_result'}[2] . ( ( defined $self->{'_no3'} ) ? '>' : '' ) )
              ),
              "\n";
            print '                     /locus_tag="', $self->{'_seqref'}->display_id, "\"\n";
            print '                     /gene="', $self->{'_result'}[5], "\"\n";
            if ( !( defined $self->{'_result'}[8] ) ) {
                print "                     /exception=\"putative or incomplete\"\n";
                my $note = $self->{'_result'}[6] . '"';
                print '                     /note="', substr( $note, 0, 51 ), "\n";
                my $i = 0;
                while ( length($note) > ( ( $i * 58 ) + 51 ) ) {
                    print '                     ', substr( $note, ( ( $i++ * 58 ) + 51 ), 58 ), "\n";
                }
            } elsif ( !( $self->{'_result'}[6] =~ /ORF1/ )
                && !( defined $self->{'_no5'} ) )
            {
                print "                     /exception=\"sequence error\"\n";
            } elsif ( defined $self->{'_minus'} ) {
                print "                     /exception=\"putatif -1 frameshift\"\n";
            } elsif ( ( defined $self->{'_no5'} ) || ( defined $self->{'_no3'} ) ) {
                print "                     /exception=\"partial sequence\"\n";
            }
            my $dna = $self->{'_result'}[7] . '"';
            print '                     /dna="', substr( $dna, 0, 52 ), "\n";
            my $i = 0;
            while ( length($dna) > ( ( $i * 58 ) + 52 ) ) {
                print '                     ', substr( $dna, ( ( $i++ * 58 ) + 52 ), 58 ), "\n";
            }
            if ( defined $self->{'_result'}[8] ) {
                print '     CDS             ',        $self->{'_result'}[0], "\n";
                print '                     /gene="', $self->{'_result'}[5], "\"\n";
                print '                     /locus_tag="', $self->{'_seqref'}->display_id, "\"\n";
                my $note = $self->{'_result'}[6] . '"';
                print '                     /note="', substr( $note, 0, 51 ), "\n";
                $i = 0;
                while ( length($note) > ( ( $i * 58 ) + 51 ) ) {
                    print '                     ', substr( $note, ( ( $i++ * 58 ) + 51 ), 58 ), "\n";
                }
                if ( ( $self->{'_result'}[6] =~ /ORF1/ ) && !( defined $self->{'_minus'} ) ) {
                    my $position3 = ( $1 + 1 )
                      if ( $self->{'_result'}[0] =~ /^\D*\d+\.\.(\d+)/ );
                    my $fs = substr(
                        $self->{'_result'}[7],
                        (
                              ( $self->{'_result'}[3] > 0 )
                            ? ( $position3 - $self->{'_result'}[1] - 15 )
                            : ( $self->{'_result'}[2] - $position3 - 15 )
                        ),
                        30
                    );
                    print '                     /inference="FS site: ', $self->_what_fs($fs), "\"\n";
                    print "                     /ribosomal_slippage\n";
                } elsif ( ( defined $self->{'_minus'} ) ) {
                    print "                     /exception=\"putatif -1 frameshift\"\n";
                } elsif ( !( defined $self->{'_no5'} ) ) {
                    print "                     /exception=\"sequence error\"\n";
                }
                print "                     /exception=\"partial sequence*\"\n"
                  if ( ( defined $self->{'_no5'} ) || ( defined $self->{'_no3'} ) );
                print "                     /codon_start=1\n"
                  unless ( defined $self->{'_noaug'} );
                print '                     /transl_table=', $self->{'_table'}, "\n";
                print '                     /product="', $self->{'_result'}[4], "\"\n";
                my $translation_issue = $self->{'_result'}[8] . '"';
                print '                     /translation="', substr( $translation_issue, 0, 44 ), "\n";
                $i = 0;

                while ( length($translation_issue) > ( ( $i * 58 ) + 44 ) ) {
                    print '                     ', substr( $translation_issue, ( ( $i++ * 58 ) + 44 ), 58 ), "\n";
                }
            }
        } elsif ( $out eq 'xml' ) {
            my $organism = shift(@args);
            print ' <sequence id="', $self->{'_seqref'}->display_id, ".seq\">\n";
            print "  <input>\n";
            print '   <organism genbank="' . $self->{'_seqref'}->display_id . "\">$organism</organism>\n"
              if ( defined $organism );
            print '   <seq type="dna" length="', length( $self->{'_seqref'}->seq ), '">', $self->{'_seqref'}->seq, "</seq>\n";
            print "  </input>\n";
            print '  <output',
              (
                ( !( $self->{'_result'}[6] =~ /ORF1/ ) && !( defined $self->{'_no5'} ) )
                ? ' comment="sequence error"'
                : ( ( defined $self->{'_minus'} ) ? ' comment="putatif -1 frameshift"' : '' )
              ),
              ' id="', $self->{'_seqref'}->display_id, "\">\n";
            print '   <gene id="', $self->{'_seqref'}->display_id, ".1\">\n";
            print '    <coord', ( defined $self->{'_no5'} ) ? ' 5prime="missing"' : '', ( defined $self->{'_no3'} ) ? ' 3prime="missing"' : '', '>',
              (
                  ( $self->{'_result'}[3] < 0 )
                ? ( 'complement(' . ( ( defined $self->{'_no3'} ) ? '<' : '' ) . $self->{'_result'}[1] . '..' . $self->{'_result'}[2] . ( ( defined $self->{'_no5'} ) ? '>' : '' ) . ')' )
                : ( ( ( defined $self->{'_no5'} ) ? '<' : '' ) . $self->{'_result'}[1] . '..' . $self->{'_result'}[2] . ( ( defined $self->{'_no3'} ) ? '>' : '' ) )
              ),
              "</coord>\n";
            print '    <name>', $self->{'_result'}[5], "</name>\n";
            print '    <seq type="dna" length="', length( $self->{'_result'}[7] ), '">', $self->{'_result'}[7], "</seq>\n";
            print "   </gene>\n";
            print '   <cds id="', $self->{'_seqref'}->display_id, ".2\">\n";
            print '    <coord', ( defined $self->{'_noaug'} ) ? ' start="unknown"' : '', ( defined $self->{'_no5'} ) ? ' 5prime="missing"' : '', ( defined $self->{'_no3'} ) ? ' 3prime="missing"' : '', '>', $self->{'_result'}[0], "</coord>\n";
            print '    <name>', $self->{'_result'}[5], "</name>\n";
            print '    <note>', $self->{'_result'}[6], "</note>\n";

            if ( defined $self->{'_result'}[8] ) {
                print '    <product>', $self->{'_result'}[4], "</product>\n";
                print '    <seq type="prt" length="', length( $self->{'_result'}[8] ), '">', $self->{'_result'}[8], "</seq>\n";
            }
            print '    <model hmm="', $self->{'_hmmoaz'}[1], '">', $self->{'_hmmoaz'}[0], "</model>\n";
            print '    <model hmm="', $self->{'_hmmfs'}[1],  '">', $self->{'_hmmfs'}[0],  "</model>\n"
              if ( $self->{'_result'}[6] =~ /ORF1/ );
            print "   </cds>\n";
            if ( ( $self->{'_result'}[6] =~ /ORF1/ ) && !( defined $self->{'_minus'} ) ) {
                my $position3 = ( $1 + 1 )
                  if ( $self->{'_result'}[0] =~ /^\D*\d+\.\.(\d+)/ );
                my $fs = substr(
                    $self->{'_result'}[7],
                    (
                          ( $self->{'_result'}[3] > 0 )
                        ? ( $position3 - $self->{'_result'}[1] - 15 )
                        : ( $self->{'_result'}[2] - $position3 - 15 )
                    ),
                    30
                );
                print '   <frameshift id="', $self->{'_seqref'}->display_id, '.3" evalue="', $self->_what_fs($fs), "\">\n";
                print '    <psite>',      substr( $fs, 12, 3 ), "</psite>\n";
                print '    <asite>',      substr( $fs, 15, 3 ), "</asite>\n";
                print '    <downstream>', substr( $fs, 18, 1 ), "</downstream>\n";
                print "   </frameshift>\n";
            }
            print "  </output>\n";
            print " </sequence>\n";
        } elsif ( $out eq 'raw' ) {    # raw format
            print $self->{'_result'}[8], "\n";
        } else {                       # fasta format
            if ( my $organism = shift(@args) ) {
                print "\n>", $organism;
            } else {
                print "\n>", $self->{'_seqref'}->display_id;
            }
            print '|', $self->{'_result'}[5], ( ( defined $self->{'_noatg'} ) ? ' NOATG' : '' ), ( ( defined $self->{'_no5'} ) ? ' NO5\'' : '' ), ( ( defined $self->{'_no3'} ) ? ' NO3\'' : '' ), ( ( $self->{'_result'}[6] =~ /ORF1/ ) ? ( ( defined $self->{'_minus'} ) ? ' PUTATIF_-1_FRAMESHIFT' : '' ) : ' SEQUENCE_ERROR' ), "\n";
            my $i = 0;
            while ( length( $self->{'_result'}[8] ) > ( $i * 80 ) ) {
                print substr( $self->{'_result'}[8], ( $i++ * 80 ), 80 ), "\n";
            }
        }
    }
}

=head2 _translation

 Title   : _translation
 Usage   : my ( $start, $end ) = $OAF_obj->_translation();
 Function: format initiation and stop codons for regex.
 Returns : array with initiation and stop codons.
 Args    : none.

=cut

sub _translation {
    my $self = shift;

    my @table        = ( 'A', 'T', 'C', 'G' );
    my @var          = ();
    my @var2         = ();
    my $var_i        = 0;
    my $var2_i       = 0;
    my $myCodonTable = Bio::Tools::CodonTable->new( -id => $self->{'_table'} );
    for ( my $i = 0 ; $i < 4 ; $i++ ) {
        for ( my $j = 0 ; $j < 4 ; $j++ ) {
            for ( my $k = 0 ; $k < 4 ; $k++ ) {
                $var[ $var_i++ ] = $table[$i] . $table[$j] . $table[$k]
                  if $myCodonTable->is_start_codon( $table[$i] . $table[$j] . $table[$k] );
                $var2[ $var2_i++ ] = $table[$i] . $table[$j] . $table[$k]
                  if $myCodonTable->is_ter_codon( $table[$i] . $table[$j] . $table[$k] );
            }
        }
    }
    @var = ('ATG') if ( $self->{'_aug'} );
    return ( join( '|', @var ), join( '|', @var2 ) );
}

=head2 get_sequence

 Title   : get_sequence
 Usage   : my @array = Bio::Tools::OAF->get_sequence( $reference );
 Function: Download sequence from Genbank safely.
 Returns : $table translation table,
           $organism organism common name,
           $seqobj PrimarySeqI objects.
 Args    : $reference (mandatory) Genbank accession number.

=cut

sub get_sequence {
    my ( $self, @args ) = @_;
    my $reference = shift(@args);

    my ( $seq, $translate, $organism );
    my $loop = 1;
    until ( ( defined $seq ) || ( $loop++ > 3 ) ) {
        my $gb = new Bio::DB::GenBank(
            -retrievaltype => 'tempfile',
            -format        => 'genbank'
        );
        eval {
            my $seqio = $gb->get_Stream_by_id($reference);
            $seq = $seqio->next_seq;
        };
    }
    if ( defined $seq->seq() ) {
        foreach my $feat ( $seq->all_SeqFeatures() ) {
            $translate = 1;
            if ( $feat->primary_tag eq 'source' ) {
                $organism = (
                    $feat->has_tag('organism')
                    ? ( $feat->each_tag_value('organism') )[0]
                    : ''
                );
            } elsif ( $feat->primary_tag eq 'CDS' ) {
                $translate = ( $feat->each_tag_value('transl_table') )[0]
                  if ( $feat->has_tag('transl_table') );
                last;
            }
        }
        return ( $translate, $organism, $seq );
    } elsif ( defined $seq ) {

        my $gb = new Bio::DB::GenBank(
            -retrievaltype => 'tempfile',
            -format        => 'fasta'
        );
        eval {
            my $seqio = $gb->get_Stream_by_id($reference);
            $seq = $seqio->next_seq;
        };
        return ( 1, 'unknown', $seq );

    }
    return undef;
}

=head2 _getBlasted

 Title   : _getBlasted
 Usage   : my @hits = $OAF_obj->_getBlasted( $bank, $evalue, $seq, $remote, $verbose );
 Function: Make a BLAST and retrieve the hits.
 Returns : Array of hits.
 Args    : $bank (optional) set the database used. If the remote NCBI
             option in set, this set the database used for the initial
             blast search. Available databases include 'refseq_rna'
             (Complete genomes and chromosomes from the NCBI Reference
             Sequence project), 'est' (Sequences from EST divisions). By
             default 'est' is used,
           $evalue (mandatory) det the E-value threshold,
           $seq (mandatory) primarySeqI object (dna or rna),
           $remote (optional) configure the program to use the remote
             BLAST server of the NCBI (1) rather than with a local
             BLAST (0). By default the local Blast is used (0),
           $verbose (optional) print more possibly useful stuff, such as
             the individual scores for each sequence.

=cut

sub _getBlasted {
    my ( $self, @args ) = @_;
    my ( $bank, $evalue, $seq, $remote, $verbose ) = @args;
    my @access;
    if ( defined $seq ) {
        if ( defined $remote && ( $remote == 1 ) ) {
            my $blasted = Bio::Tools::BlastTools->new_query( $bank, 'tblastn', $evalue, 10000 );
            if ( $blasted->submit( $seq, undef, ( $evalue < 1e-10 ) ? 7 : undef ) ) {
                my $hits = Bio::Tools::BlastTools->new_hits(1);
                print STDERR 'Remote Blast (', $blasted->rid, ")\n" if ($verbose);
                if ( $hits->retrieve($blasted) ) {
                    for ( my $i = 0 ; $i < $hits->hits ; $i++ ) {
                        push @access, ( split( /\./, ( split( /\|/, $hits->name($i) ) )[3] ) )[0] . '|' . $hits->hitstart($i) . '|' . $hits->hitend($i) . '|' . $hits->hitstrand($i) . '|' . $hits->hitframe($i);
                    }
                }
            }
        } else {
            print STDERR "Local Blast\n" if ($verbose);
            my $blasted = Bio::Tools::BlastLocalTools->new_query( $bank, 'tblastn', $evalue, 10000 );
            if ( $blasted->submit($seq) ) {
                for ( my $i = 0 ; $i < $blasted->hits ; $i++ ) {
                    push @access, ( split( /\./, ( split( /\|/, $blasted->name($i) ) )[1] ) )[0] . '|' . $blasted->hitstart($i) . '|' . $blasted->hitend($i) . '|' . $blasted->hitstrand($i) . '|' . $blasted->hitframe($i);
                }
            }
        }
    }
    return @access;
}

=head2 ortholog

 Title   : ortholog
 Usage   : my $bool = Bio::Tools::OAF->ortholog( $file, $evalue, $format, $aug, $remote, $hmm, $fs, $bank, $filter, $verbose );
 Function: Full BLAST search...
 Returns : A boolean.
 Args    : $file (mandatory) name of the file of protein sequence(s) (FASTA file format),
           $evalue (optional) set the E-value threshold,
           $format (optional) output format: 'fasta', 'genbank', 'raw' or 'xml' style,
           $aug (optional) use other start codon than AUG (default 0),
           $remote (optional) configure the program to use the remote
             BLAST server of the NCBI (1) rather than with a local
             BLAST (0). By default the local Blast is used (0),
           $hmm (optional) path to hmm profiles by default OAF looks at
             ./oaz.hmm,
           $fs (optional) path to frameshift hmm profiles by default OAF
             looks at ./oaz.fs.hmm,
           $bank (optional) set the database used. If the remote NCBI
             option in set, this set the database used for the initial
             blast search. Available databases include 'refseq_rna'
             (Complete genomes and chromosomes from the NCBI Reference
             Sequence project), 'est' (Sequences from EST divisions). By
             default 'est' is used,
           $filter (optional) set the filter. Normally OAF.pl show all
             sequences identified whether the sequence may be the
             putative, because of a truncated EST, a sequence error or a
             missing -1 frameshift. By default there is no filter (0), you
             can filter for the non-ambiguous hits only (1) or the
             ambiguous hits only (2),
           $verbose (optional) print more possibly useful stuff, such as
             the individual scores for each sequence.

=cut

sub ortholog {
    my ( $self, @args ) = @_;
    my ( $file, $evalue, $format, $aug, $remote, $hmm, $fs, $bank, $filter, $verbose ) = @args;
    $remote  = 0       unless ( defined $remote );
    $filter  = 0       unless ( defined $filter );
    $verbose = 0       unless ( defined $verbose );
    $format  = 'fasta' unless ( defined $format );
    $bank    = 'est'   unless ( defined $bank );
    my @access;

    my $inseq = Bio::SeqIO->new( '-file' => "<$file", '-format' => 'fasta' );
    while ( my $seq = $inseq->next_seq ) {
        print ">", $seq->display_id if ($verbose);
        @access = keys %{ { map { $_ => 1 } ( @access, $self->_getBlasted( $bank, 1e-5, $seq, $remote, $verbose ) ) } };
    }
    if ( $#access >= 0 ) {
        @access = sort @access;
        my @accession;
        for ( my $i = 0 ; $i <= $#access ; $i++ ) {
            my ( $chrom, $hitstart, $hitend, $hitstrand, $hitframe ) =
              split( /\|/, $access[$i] );
            my $gloops = 0;
            do {
                $gloops = 0;
                my ( $chrom_n, $hitstart_n, $hitend_n, $hitstrand_n, $hitframe_n ) = split( /\|/, $access[ $i + 1 ] )
                  if ( defined $access[ $i + 1 ] );

                #frameshift or echo ?
                if (
                       ( defined $chrom_n )
                    && ( $chrom eq $chrom_n )
                    && (   ( abs( $hitstart_n - $hitend ) < 1500 )
                        || ( abs( $hitend_n - $hitend ) < 1500 )
                        || abs( $hitstart_n - $hitend ) < 1500 )
                    && ( $hitstrand == $hitstrand_n )
                  )
                {
                    $i++;
                    $hitend   = $hitend_n   if ( $hitend < $hitend_n );
                    $hitstart = $hitstart_n if ( $hitstart > $hitstart_n );
                    $gloops   = 1;
                }
            } until ( $gloops == 0 );
            push @accession, ( $chrom . '|' . $hitstart . '|' . $hitend . '|' . $hitstrand . '|' . $hitframe );
        }
        @access = @accession;
        undef @accession;
    }

    print STDERR "\nRetrieve hits (", ( $#access + 1 ), ")\n" if ($verbose);
    if ( $#access >= 0 ) {
        my ( $translate, $organism ) = ( 0, '' );
        my ( $seq_chr, $last, @cds_features );
        @access = sort @access;
        for ( my $i = 0 ; $i <= $#access ; $i++ ) {
            my ( $chrom, $hitstart, $hitend, $hitstrand, $hitframe ) =
              split( /\|/, $access[$i] );
            unless ( defined $last && $last eq $chrom ) {
                ( $translate, $organism, $seq_chr ) = $self->get_sequence($chrom);
                @cds_features =
                  grep { $_->primary_tag eq 'CDS' } $seq_chr->get_SeqFeatures;
                print STDERR "\r> Load ", $seq_chr->accession_number, " ($organism)\n  "
                  if ($verbose);
                $last = $chrom;
            }
            next unless ( $seq_chr->molecule =~ /RNA/ );
            my $seqstart = $hitstart - 350;
            $seqstart = 1 if ( $seqstart < 1 );
            my $seqend = $hitend + 350;
            $seqend = $seq_chr->length() if ( $seqend > $seq_chr->length() );
            my $seq = Bio::Seq->new(
                -seq      => $seq_chr->subseq( $seqstart, $seqend ),
                -alphabet => 'dna',
                -id       => $chrom
            );
            my $OAF_obj = $self->new( $seq, $translate, $aug, $hmm, $fs );
            $OAF_obj->{'_verbose'} = 1 if ($verbose);
            $OAF_obj->{'_filter'} = $filter;

            if ( $OAF_obj->find( $evalue, $hitstrand, $seqstart, $seqend ) ) {
                $OAF_obj->show($format);
            } elsif ( $verbose && ( defined $OAF_obj->{'_hmmoaz'}[0] ) && ( $OAF_obj->{'_hmmoaz'}[0] < $OAF_obj->{'_evalue'} ) ) {
                $OAF_obj->show($format);
            }
        }
        return 1;
    }
    return 0;
}

# and that's all the module
1;
