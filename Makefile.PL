# $Id: Makefile.pm,v 1.0.6 2008/01/25 02:53:11 Europe/Dublin $

require 5.005;
use ExtUtils::MakeMaker ;
use File::Path 'rmtree';
use IO::File;
use Config;

$NAME     = 'Bio';
$DISTNAME = 'oaf';
$VERSION  = '1.0.6';

WriteMakefile(
   VERSION => $VERSION,
   NAME	=> $NAME,
   DISTNAME    => $DISTNAME,
   'dist'      => {  COMPRESS     => 'bzip2',
                     SUFFIX       => '.bz2',
                     DIST_DEFAULT => 'all tardist',
                  },
   ($] ge '5.005')  ? (
    'AUTHOR'   => 'Michael Bekaert (mbekaert@gmail.com)',
    'ABSTRACT' => 'Ornitine decarboxylase Antizyme Finder (OAF)',
   ) : (),
   DIR         => [  ],
);
