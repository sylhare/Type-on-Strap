---
published: true
tags:
  - networking
layout: post
---
# Layer 3 THE NETWORK LAYER

unit : datagram

### Virtual circuit network: 
it represented for incomming/outgoing interface in logic on each port of source

example one packet can come from multiple interface on a router.After that the router send out another interface with other VC number

+![virtual circuit table]({{site.baseurl}}/assets/img/virtual_circuit_table.png)

## Inside a Router 

+![router_architecture]({{site.baseurl}}/assets/img/router_architecture.png)

### Input port processing
in this place. Router decapsultation, lookup/forwarding  or queing to switch fabric ( prepare next step)

### switch fabric
in his point . there are 3 way to Router can switch message from input to output 

+![switching techniques]({{site.baseurl}}/assets/img/switching_techniques.png)
1. via memory . it is controled by processor . speed quite slow from milisecond -> hundred milisecond
2. via bus . hardware with configuration by software before . quicker than type one but bus can service 1 message at a time
3. crossbar . conenct input and output with many intersection => have more way one message can go from SRC to DES quickest this technique support send parallel multiple message in case all of them differ DES

### Output processing 

At this point. Router selecting and de-queueing packet for tranmission and performing the needed link-layer and physical-layer transmission functions
