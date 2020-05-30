## A New Post
The transport layer at the server notes the following four values in the con-nection-request segment: 
(1) the source port number in the segment, 
(2) the IPaddress of the source host, 
(3) the destination port number in the segment
(4) its own IP address. 
The newly created connection socket is identified by thesefour values; all subsequently arriving segments whose source port, source IPaddress, destination port, and destination IP address match these four values willbe demultiplexed to this socket. With the TCP connection now in place, the clientand server can now send data to each other.

The server host may support many simultaneous TCP connection sockets, witheach socket attached to a process, and with each socket identified by its own four-tuple. When a TCP segment arrives at the host, all four fields (source IP address,source port, destination IP address, destination port) are used to direct (demultiplex)the segment to the appropriate socket

For TCP, nmap sequentiallyscans ports, looking for ports that are accepting TCP connections. For UDP, nmapagain sequentially scans ports, looking for UDP ports that respond to transmitted UDP segments. In both cases, nmap returns a list of open, closed, or unreachableports. A  host running nmap can attempt to scan any target host anywherein theInternet

UDP require a minimum seding rate & without hand shaking step ( without the delay)
UDP :No connection establishment
UDP :No Connection state 

UDP, on the other hand, doesnot maintain connection state and does not track any of these parameters. For thisreason, a server devoted to a particular application can typically support manymore active clients when the application runs over UDP rather than TCP.•Small packet header overhead.

The TCP segment has 20 bytes of header over-head in every segment, whereas UDP has only 8 bytes of overhead
