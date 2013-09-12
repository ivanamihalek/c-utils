#! /sur/bin/perl -w

while (<>){
    chomp;
    push @val , split '';
} 
foreach $i ( 0 .. $#val) {
    printf (" %3d %3d %3d \n", int ($i/50), $i%50,  $val[$i]);
}
