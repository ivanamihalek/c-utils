#! /bin/bash


# don't foget: 
# valgrind --leak-check=yes ./daddy /databases/db.txt

# newline in the input buffer
echo
echo "test: newline in the input buffer"
gcc -o pinger  pinger.c  -DNWLN
./pinger "neighbors 1  57572  "

# not waiting for response
echo "test: not waiting for response"
gcc -o pinger  pinger.c  -DSKIP_RECV
./pinger "neigbors 1  57572  "

# no method name
echo "test: missing method name"
gcc -o pinger  pinger.c 
./pinger " 57572 7049"

echo "test: bad method name"
# broken method name
./pinger "bnmetq 57572 7049"

echo "test: space before method name"
# space before  method name
./pinger " neighbors 1  57572"


# order out of range
echo "test: order out of range"
./pinger " neighbors -1  57572 "
./pinger " neighbors 10  57572"
./pinger " neighbors 1y  57572"

# input that is not number(s)
echo
echo "test: input that is not number(s) "
./pinger " neighbors 1  57572abg "

# very short input 
echo
echo "test: very short input "
./pinger "neighbors 1   "
./pinger "path  1  57572 "

# very long  input
echo
echo "test: very long input"
gcc -o pinger  pinger.c -DLONG_INPUT
./pinger "neighbors 1 "

# very too long  input
echo
echo "test: too long input"
gcc -o pinger  pinger.c -DTOO -DLONG_INPUT
./pinger "neighbors 1 "


# newline in the input buffer
echo
echo "test: shortest path"
gcc -o pinger  pinger.c  
./pinger "neighbors 0  57572  4079 "
