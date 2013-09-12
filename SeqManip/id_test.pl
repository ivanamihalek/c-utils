#! /usr/gnu/bin/perl
use SeqManip;




$file1 = "1dqwA.input";

open ( FILE1, "<$file1") 
|| die "open fail $file1"; 
$ok = "n";
$length = 30;
$percent_id = 15.0;
while ( <FILE1>) {

    next if />/;
    print $_;
    chomp;
    $retval = SeqManip::identity ( $percent_id, $_, $length, $ok);
    print "ok: $ok    retval: $retval\n";
    print "\n\n"; 
}
