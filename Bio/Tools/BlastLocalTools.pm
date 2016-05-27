## $Id: BlastLocalTools.pm,v 1.10 2007/03/17 16:56:27 Europe/Dublin $
#
# BioPerl module for Bio::Tools::BlastLocalTools
#
# 'Blast made easy' shared package
# Copyright 2005-2007 Bekaert M <mbekaert@gmail.com>
#
# This work is licensed under the Creative Commons Attribution-
# Noncommercial-Share Alike 3.0 License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BlastLocalTools - Local Blast made easy

=head1 SYNOPSIS

Take a sequence object from, say, an inputstream, and return
an object including all BLAST hits.

Creating the BlastTools object, eg:

  my $request=Bio::Tools::BlastLocalTools->new_query();

submit a request and obtain results, eg:

  my $inputstream = Bio::SeqIO->new(-file => 'seqfile', -format => 'fasta');
  my $seqobj = $inputstream->next_seq();
  $request->submit($seqobj);
  for(my $i=0; $i< $request->hits; $i++) {
    print $request->name($i);
  }


 =head1 DESCRIPTION

Bio::Tools::BlastLocalTools is a featherweight object for an easy
access to StandAlone BLAST. Returned object include location,
sequence and statistic for each hit found.

See Synopsis above for object creation code.

=head1 DRIVER SCRIPT

 #!/usr/bin/perl -w
 use strict;

 use Bio::Seq;
 use Bio::Tools::BlastLocalTools;

 my $inseq = Bio::SeqIO->new('-file' => "< $yourfile", -format => 'Fasta');
 while (my $seq = $inseq->next_seq) {
   my $request=Bio::Tools::BlastLocalTools->new_query();
   if ($request->submit($seqobj)) {
     print '> ',$request->hits," hit(s)\n";
   } else {
     print "  no hit!\n";
   }
 }

=head1 REQUIREMENTS

To use this module you need Bioperl (L<http://www.bioperl.org/>) modules and
a local Blast version (L<http://www.ncbi.nlm.nih.gov/BLAST/download.shtml>).

=head1 FEEDBACK

User feedback is an integral part of the evolution of this
modules. Send your comments and suggestions preferably
to author.

=head1 AUTHOR

B<Michael Bekaert> (mbekaert@gmail.com)

Address:
     School of Biology & Environmental Sciences
     University College Dublin
     Belfield, Dublin 4
     Dublin
     Ireland

=head1 SEE ALSO

perl(1), blastall(1) and bioperl web site

=head1 LICENSE

Copyright 2005-2007 - Michael Bekaert

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::BlastLocalTools;

use vars qw(@ISA $VERSION);
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;

use Bio::Tools::Run::StandAloneBlast;
use Bio::Search::Hit::BlastHit;

$VERSION = '1.10';

@ISA = qw(Bio::Root::Root Bio::Root::IO);

=head2 new_query

 Title   : new_query
 Usage   : my $query = Bio::Tools::BlastLocalTools->new_query($dbname,$prog,$evalue,$max);
 Function: Initialize Query object.
 Returns : Query object.
 Args    : $dbname (optional) database used (default nr),
           $prog (optional) BLAST programme used (default blastn),
           $evalue (optional) Set the E-value threshold (default 1e-100),
           $max (optional) Maximum sequence wait (default 100).

=cut

sub new_query {
    my ( $self, @args ) = @_;
    my ( $dbname, $prog, $evalue, $max ) = @args;
    $self = {};
    $self->{'_dbname'} = ( ( defined $dbname ) ? $dbname : 'est' );
    $self->{'_prog'}   = ( ( defined $prog )   ? $prog   : 'blastn' );
    $self->{'_zvalue'} = ( ( defined $evalue ) ? $evalue : '1e-100' );
    $self->{'_max'}    = ( ( defined $max )    ? $max    : '100' );
    $self->{'_seq'}   = undef;
    $self->{'_table'} = undef;
    $self->{'_hits'}  = 0;
    bless($self);
    return $self;
}

=head2 hits

 Title   : hits
 Usage   : my $hits = $hits->hits();
 Function: Return the number of results/hits.
 Returns : $hits Integer.
 Args    : none.

=cut

sub hits {
    my $self = shift;
    return $self->{'_hits'};
}

=head2 name

 Title   : name
 Usage   : my $r = $hits->name(1);
 Function: Get the name of the first (1) result.
 Returns : Return a string
 Args    : Index of the result.

=cut

sub name {
    my $self = shift;
    if (@_) { return $self->{'_name'}[shift]; }
}

=head2 description

 Title   : description
 Usage   : my $r = $hits->description(1);
 Function: Get the description of the first (1) result.
 Returns : Return a string
 Args    : Index of the result.

=cut

sub description {
    my $self = shift;
    if (@_) { return $self->{'_description'}[shift]; }
}

=head2 evalue

 Title   : evalue
 Usage   : my $r = $hits->evalue(1);
 Function: Get the e-value of the first (1) result.
 Returns : Return a float
 Args    : Index of the result.

=cut

sub evalue {
    my $self = shift;
    if (@_) { return $self->{'_evalue'}[shift]; }
}

=head2 bits

 Title   : bits
 Usage   : my $r = $hits->bits(1);
 Function: Get the number of bits of the first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub bits {
    my $self = shift;
    if (@_) { return $self->{'_bits'}[shift]; }
}

=head2 qlength

 Title   : qlength
 Usage   : my $r = $hits->qlength(1);
 Function: Get the length of the query sequence used for the
            first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub qlength {
    my $self = shift;
    if (@_) { return $self->{'_qlength'}[shift]; }
}

=head2 querystart

 Title   : querystart
 Usage   : my $r = $hits->querystart(1);
 Function: Get the 5' position of the query sequence used for the
           first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub querystart {
    my $self = shift;
    if (@_) { return $self->{'_querystart'}[shift]; }
}

=head2 queryend

 Title   : queryend
 Usage   : my $r = $hits->queryend(1);
 Function: Get the  3' position of the query sequence used for the
           first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub queryend {
    my $self = shift;
    if (@_) { return $self->{'_queryend'}[shift]; }
}

=head2 querygaps

 Title   : querygaps
 Usage   : my $r = $hits->querygaps(1);
 Function: Get the gap number of the query sequence for the
           first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub querygaps {
    my $self = shift;
    if (@_) { return $self->{'_querygaps'}[shift]; }
}

=head2 querystrand

 Title   : querystand
 Usage   : my $r = $hits->querystrand(1);
 Function: Get the strand of the query sequence for the
           first (1) result (-1 Crick, 1 Watson).
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub querystrand {
    my $self = shift;
    if (@_) { return $self->{'_querystrand'}[shift]; }
}

=head2 queryframe

 Title   : queryframe
 Usage   : my $r = $hits->queryframe(1);
 Function: Get the frame of the query sequecence
           for the first (1) result (proteins only).
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub queryframe {
    my $self = shift;
    if (@_) { return $self->{'_queryframe'}[shift]; }
}

=head2 hitprimary

 Title   : hitprimary
 Usage   : my $r = $hits->hitprimary(1);
 Function: Get the ? of the first (1) result.
 Returns : Return a string
 Args    : Index of the result.

=cut

sub hitprimary {
    my $self = shift;
    if (@_) { return $self->{'_hitprimary'}[shift]; }
}

=head2 hitsource

 Title   : hitsource
 Usage   : my $r = $hits->hitsource(1);
 Function: Get the ? of the first (1) result.
 Returns : Return a string
 Args    : Index of the result.

=cut

sub hitsource {
    my $self = shift;
    if (@_) { return $self->{'_hitsource'}[shift]; }
}

=head2 hitstart

 Title   : hitstart
 Usage   : my $r = $hits->hitstart(1);
 Function: Get the 5' position of the first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub hitstart {
    my $self = shift;
    if (@_) { return $self->{'_hitstart'}[shift]; }
}

=head2 hitend

 Title   : hitend
 Usage   : my $r = $hits->hitend(1);
 Function: Get the 3' position of the first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub hitend {
    my $self = shift;
    if (@_) { return $self->{'_hitend'}[shift]; }
}

=head2 hitgaps

 Title   : hitgaps
 Usage   : my $r = $hits->hitgaps(1);
 Function: Get the gap number of the first (1) result.
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub hitgaps {
    my $self = shift;
    if (@_) { return $self->{'_hitgaps'}[shift]; }
}

=head2 hitstrand

 Title   : hitstrand
 Usage   : my $r = $hits->hitstrand(1);
 Function: Get the strand of the first (1) result
           (-1 Crick, 1 Watson).
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub hitstrand {
    my $self = shift;
    if (@_) { return $self->{'_hitstrand'}[shift]; }
}

=head2 hitframe

 Title   : hitframe
 Usage   : my $r = $hits->hitframe(1);
 Function: Get the frame of the first (1) result (proteins only).
 Returns : Return an integer
 Args    : Index of the result.

=cut

sub hitframe {
    my $self = shift;
    if (@_) { return $self->{'_hitframe'}[shift]; }
}

=head2 hitseq

 Title   : hitseq
 Usage   : my $r = $hits->hitseq(1);
 Function: Get the sequence of the first (1) result.
 Returns : Return a string
 Args    : Index of the result.

=cut

sub hitseq {
    my $self = shift;
    if (@_) { return $self->{'_hitseq'}[shift]; }
}

=head2 submit

 Title   : submit
 Usage   : my $statut = $query->submit( $seqobj, $table, $word );
 Function: submit a Blast search.
 Returns : submission status.
 Args    : $seqobj (mandatory) PrimarySeqI object,
           $table (optional) translation table number (default 1),
           $word (optional) word size.

=cut

sub submit {
    my ( $self, @args ) = @_;
    my ( $seqobj, $table, $word ) = @args;

    if ( ( defined $seqobj ) && ( defined $seqobj->seq ) ) {
        $self->{'_seq'} = $seqobj;
        $self->{'_table'} = ( ( defined $table ) ? $table : 1 );
        my @params = (
            -program  => $self->{'_prog'},
            -datadase => $self->{'_dbname'},
            -e        => $self->{'_zvalue'}
        );
        my $factory      = Bio::Tools::Run::StandAloneBlast->new(@params);
        my $blast_report = $factory->blastall( $self->{'_seq'} );
        if ( ref($blast_report) ) {
            my $result = $blast_report->next_result;
            $self->{'_hits'} = 0;
            while ( my $hit = $result->next_hit ) {
                while (( my $hsp = $hit->next_hsp )
                    && ( $self->{'_hits'} <= $self->{'_max'} ) )
                {
                    $self->{'_name'}[ $self->{'_hits'} ]        = $hit->name;
                    $self->{'_description'}[ $self->{'_hits'} ] = $hit->description();
                    $self->{'_evalue'}[ $self->{'_hits'} ]      = $hsp->evalue();
                    $self->{'_bits'}[ $self->{'_hits'} ]        = $hsp->bits;
                    $self->{'_qlength'}[ $self->{'_hits'} ]     = $hsp->length;
                    $self->{'_querystart'}[ $self->{'_hits'} ]  = $hsp->query->start;
                    $self->{'_queryend'}[ $self->{'_hits'} ]    = $hsp->query->end;
                    $self->{'_hitprimary'}[ $self->{'_hits'} ]  = $hsp->hit->primary_tag;
                    $self->{'_hitsource'}[ $self->{'_hits'} ]   = $hsp->hit->source_tag;
                    $self->{'_hitstart'}[ $self->{'_hits'} ]    = $hsp->hit->start;
                    $self->{'_hitend'}[ $self->{'_hits'} ]      = $hsp->hit->end;
                    $self->{'_querygaps'}[ $self->{'_hits'} ]   = $hsp->gaps('query');
                    $self->{'_hitgaps'}[ $self->{'_hits'} ]     = $hsp->gaps('sbjct');
                    $self->{'_querystrand'}[ $self->{'_hits'} ] = $hsp->strand('query');
                    $self->{'_hitstrand'}[ $self->{'_hits'} ]   = $hsp->strand('sbjct');
                    $self->{'_hitframe'}[ $self->{'_hits'} ]    = $hsp->hit->frame;
                    $self->{'_queryframe'}[ $self->{'_hits'} ]  = $hsp->query->frame;
                    $self->{'_hitseq'}[ $self->{'_hits'} ]      = $hsp->hit_string;
                    $self->{'_hits'}++;
                }
            }
            return 1;
        }
    }
}

# and that's all the module
1;
