#! /bin/bash
# newline in the input buffer
echo "test: newline in the input buffer"
gcc -o pinger  daddy_pinger.c  -lpthread  -ligraph    -I/usr/local/include/igraph -DNWLN
./pinger "neighbors 1  57572  "

# not waiting for response
echo "test: not waiting for response"
gcc -o pinger  daddy_pinger.c  -lpthread  -ligraph    -I/usr/local/include/igraph -DSKIP_RECV
./pinger "neigbors 1  57572  "

# no method name
echo "test: missing method name"
gcc -o pinger  daddy_pinger.c  -lpthread  -ligraph    -I/usr/local/include/igraph
./pinger " 57572 7049"

echo "test: bad method name"
# broken method name
./pinger "bnmetq 57572 7049"

echo "test: space before method name"
# space before  method name
./pinger " neighbors 1  57572"

# very short input 
echo "test: very short input "
./pinger "neighbors 1   "
./pinger "path  1  57572 "

# input that is not number(s)
echo "test: input that is not number(s) "
./pinger " neighbors 1  57572abg "

# very long  input

