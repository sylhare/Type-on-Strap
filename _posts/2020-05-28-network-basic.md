---
published: false
---
## Theory , 4 tuples value
The transport layer at the server notes the following four values in the con-nection-request segment: 
(1) the source port number in the segment, 
(2) the IPaddress of the source host, 
(3) the destination port number in the segment
(4) its own IP address. 
The newly created connection socket is identified by thesefour values; all subsequently arriving segments whose source port, source IPaddress, destination port, and destination IP address match these four values willbe demultiplexed to this socket.

When a TCP segment arrives at the host, all four fields (source IP address,source port, destination IP address, destination port) are used to direct (demultiplex)the segment to the appropriate socket

## scanning port
For TCP, nmap sequentiallyscans ports, looking for ports that are accepting TCP connections. For UDP, nmapagain sequentially scans ports, looking for UDP ports that respond to transmitted UDP segments. In both cases, nmap returns a list of open, closed, or unreachableports. A host running nmap can attempt to scan any target host anywherein theInternet

## difference between TCP & UDP
UDP require a minimum seding rate & without hand shaking step ( without the delay)
UDP :No connection establishment
UDP :No Connection state 

UDP, on the other hand, doesnot maintain connection state and does not track any of these parameters. For this reason, a server devoted to a particular application can typically support manymore active clients when the application runs over UDP rather than TCP. Small packet header overhead.

The TCP segment has 20 bytes of header over-head in every segment, whereas UDP has only 8 bytes of overhead

Blocking UDP traffic for security reasons, TCP becomes an increasingly attractive protocol for streaming media transport.

## principles of reliable data transfer 
reliable is implement on diffence layers. each layer has own implementation because may be below it don't have any implementaion reliable.

### ARQ Automaic Repeat reQuest protocols
on receiver side:

> error detection ( via checksum field)

> receiver feedback (ACK = 1 or NAK = 0) 

Retransimition : A packet that is received in error at the receiver will be retrans-mitted by the sender
TCP have sequence number to indicate piece of packet number

# conclusion 
 Check-sums, sequence numbers, timers, and positive and negative acknowledgment pack-ets each play a crucial and necessary role in the operation of the protocol
 
 ### impovement in tcp 
Go-Back-N (GBN) protocol :  the sender is allowed to transmit multiple packets(when available) without waiting for an acknowledgment,
Selective Repeat (SR)    : make sender only resend with specific lost packet instead of all packet in a window size

![summary of reliable data transfer](summary_reliable_data_transfer.png)
