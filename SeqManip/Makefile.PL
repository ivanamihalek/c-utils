#! /usr.gnu/bin/perl

use ExtUtils::MakeMaker;

WriteMakefile(
'NAME'     => 'SeqManip',
'OBJECT'   => 'identity.o  SeqManip_wrap.o'
			);
