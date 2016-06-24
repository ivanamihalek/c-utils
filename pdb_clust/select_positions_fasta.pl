#! /usr/bin/perl -w
# pdb clustering prog expects the input in the
# format <aa_short>\s*<position> one per line, no space
# this script will produce it from the fasta file (sequence)
# an a list of positions, one per line

@ARGV >= 2  || 
    die "usage: $0  <fasta file>  <positions file>\n";

($fasta, $pos_file) = @ARGV;
(-e $fasta) || die "$fasta not found\n";
(-e $pos_file) || die "$pos_file not found\n";

$seq = `grep -v '>' $fasta`;
$seq =~ s/\n//g;


@pos = split '\n', `cat $pos_file`;
for $pos (@pos) {
    $val = substr $seq, $pos-1, 1;
    print "$val   $pos\n";
}
