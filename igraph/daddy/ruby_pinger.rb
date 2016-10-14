#!/usr/bin/env ruby
require 'socket'

socket = UNIXSocket.new("/tmp/igraph_daddy")

socket.write "neighbors  1 57572  7049 2263 "


data = ""
recv_length = 1024
while  tmp = socket.recv(recv_length)
    data += tmp
    break if tmp.length < recv_length
end
puts data.rstrip

socket.close
