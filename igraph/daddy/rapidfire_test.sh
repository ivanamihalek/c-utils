#! /bin/bash

cdc42=2744
tgfbr3=16386
dock6=4520
rs1=14036

echo; echo "==================================="
echo "test: bombarding the socket with requests"
gcc -o pinger  pinger.c
for  i in {1..50}
do
    echo "request $i"
    ./pinger " neighbors 10  $cdc42 "  & 
    #../pinger " neighbors 10  $cdc42 " |  cut -c-20 
    #./pinger "path  0  $cdc42 $tgfbr3  "  
    #./pinger "path  0  $cdc42 $tgfbr3  "  | tee test.out
    #grep  -q 'timed out' test.out &&  exit 
    
done	  






