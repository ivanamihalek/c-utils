#! /bin/bash

cdc42=2744
tgfbr3=16386
dock6=4520

# don't foget: 
# valgrind --leak-check=yes ./daddy /databases/db.txt

# newline in the input buffer
echo
echo "test: newline in the input buffer"
gcc -o pinger  pinger.c  -DNWLN
./pinger "neighbors 1   $cdc42 "

# not waiting for response
echo "test: not waiting for response"
gcc -o pinger  pinger.c  -DSKIP_RECV
./pinger "neigbors 1  $cdc42 " 

# no method name
echo "test: missing method name"
gcc -o pinger  pinger.c 
./pinger  "$cdc42 $tgfbr3"

echo "test: bad method name"
# broken method name
./pinger "bnmetq   $cdc42 $tgfbr3 "

echo "test: space before method name"
# space before  method name
./pinger " neighbors 1  $cdc42 $tgfbr3 "


# order out of range
echo "test: order out of range"
./pinger " neighbors -1   $cdc42 "
./pinger " neighbors 10  $cdc42 "
./pinger " neighbors 1y  $cdc42 " 

# input that is not number(s)
echo
echo "test: input that is not number(s) "
./pinger " neighbors 1  57572abg "

# very short input 
echo
echo "test: very short input "
./pinger "neighbors 1   "
./pinger "path  1  $cdc42"

# very long  input
echo
echo "test: very long input"
gcc -o pinger  pinger.c -DLONG_INPUT
./pinger "neighbors 1 "

#  too long  input
echo
echo "test: too long input"
gcc -o pinger  pinger.c -DTOO -DLONG_INPUT
./pinger "neighbors 1 "


# shortest path(s)
echo
echo "test: shortest path"
gcc -o pinger  pinger.c  
./pinger "path  0  $cdc42 $tgfbr3  "  

# node without neighbors
echo
echo "test: node with no neighbors - zero hops away"
gcc -o pinger  pinger.c  
./pinger "neighbors  0   $dock6"
echo
echo "test: node with no neighbors - zero hops away"
gcc -o pinger  pinger.c  
./pinger "neighbors  3  $dock6" 




