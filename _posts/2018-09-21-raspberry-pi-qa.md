---
layout: post
title: Raspberry Pi for QA
description: The design and implementation of multiple QA services (Network tweaking, VPN, Packet analysis,..) on a NodeJS server on a Raspberry Pi
author-id: "galera"
categories: [raspberry-pi]
tags: [linux,nodejs,raspberry-pi,qa,testing]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/rpiqa/featured-image.jpg"
thumbnail: "assets/img/posts/rpiqa/featured-image.jpg"
image: "assets/img/posts/rpiqa/featured-image.jpg"
redirect_from:
    - /2018/09/21/raspberry-pi-for-qa/
---
<p>On a previous company the QA team needed to perform some scenarios that were difficult to reproduce, for instance force a buffering event. Their setup require deep technical understanding of Linux. In order to ease its world, I decided to create a neat service running on a Raspberry Pi.</p>
<p><!--more--></p>
<h2>The problems:</h2>
<p>Concretely, the scenarios that the QA needed to cover were:</p>
<h3><a href="#network-performance"><em>Restrict the network performance</em></a></h3>
<p>The test that the QA team were performing, included provoking buffering issues in media players. They found issues on how this could be achieved: in Chrome it could be achieved through Developer Tools, but what about mobile devices, or even worst: some weird embedded media players?</p>
<h3><a href="#geoblock"><em>Access Geo-blocked content</em></a></h3>
<p>Since the company had customers world-wide and some of them used Geo-blocked content, they need to access those contents through the use of VPNs. This required spending some time configuration the VPN clients on their side</p>
<h3><a href="#ipblock"><em>IP Blocking</em></a></h3>
<p>There were some scenarios (I cannot remember right now) that required to block the connection to some IP. This is really easy to do on a UNIX machine with iptables, but good luck doing that on Windows.</p>
<h3><a href="#dnsspoof"><em>Modify DNS records</em></a></h3>
<p>It were some scenario where they needed to change the DNS records. For instance: www.google.com -&gt; 192.168.1.100. I don't remember the rationale to this requirement :(</p>
<h3><a href="#httpssniffer"><i>Analyse</i><em> HTTPS traffic</em></a></h3>
<p>In the majority of the environments the connections with the company server's were made with HTTPS. This caused a little bit of headache while analysing the HTTP traffic.</p>
<h2>The solutions ...</h2>
<p>Since most of the scenarios require some networking tweaks, the obvious decision was Linux, even better: Raspberry Pi. The project consist in a NodeJS Express application that executed scripts on the RaspberryPi.</p>
<p>The raspberry PI have two network interfaces: the LAN and the WiFi. The configuration is setup in the interfaces file:</p>

```
auto lo
iface lo inet loopback
auto eth0
auto wlan0
#static ethernet conf
iface eth0 inet static
address 192.168.1.100
netmask 255.255.255.0
gateway	192.168.1.1
dns-nameservers 8.8.4.4
iface wlan0 inet static
address 192.168.150.1
netmask 255.255.255.0
```

<p>In the WiFi interface hostapd service is configured in order to serve a WiFi connection. This will configure the traffic outgoing traffic from wlan0 to eth0.</p>
<h3>Playing with tc and ifb</h3>
<p>The network performance restrictions can be applied using the tc command and the ifb module on Linux. This module redirects the traffic from one real network interface to a virtual one. When the traffic passes through the ifb0 interface the token bucket (<a href="https://en.wikipedia.org/wiki/Token_bucket#Hierarchical_token_bucket">htb</a>) applies the network configuration</p>

```
#Enable ifb module and setup the ifb0 interface
modprobe ifb numifbs=1 &amp;&amp; ip link set dev ifb0 up
#Create a device that redirects all the traffic
#from eth0 to ifb0
tc qdisc add dev eth0 handle ffff: ingress &amp;&amp; \
tc filter add dev eth0 parent ffff: protocol ip u32 \
match u32 0 0 action mirred egress redirect dev ifb0
#Modify the token bucket configuration
tc qdisc add dev ifb0 root handle 1: htb default 10 &amp;&amp; \
tc class add dev ifb0 parent 1: classid 1:1 htb rate 1mb &amp;&amp; \
tc class add dev ifb0 parent 1:1 classid 1:10 htb rate 1mb
```

<h3>OpenVPN and Iptables</h3>
<p>In order to be able to avoid the Geo-blocking of contents, the company provide us a commercial VPN. This consisted on a series of OpenVPN configuration scripts for multiple countries. Therefore, the app only needs to call openvpn:</p>

```
openvpn ALBANIA-TCP.ovpn
```

<p>It's really easy to block an IP on a Linux box, the app only needs to call the iptables, for example:</p>

```
iptables -I FORWARD -s 192.168.150.0/24 -d 8.8.8.8  -j DROP
```

<p>This snippet blocks the outgoing connections to 8.8.8.8 that goes out from the WiFi network</p>
<h3>Networking stuff ...</h3>
<p id="dnsspoof">Regarding DNS Spoofing, it's a little bit trickier, however the <strong>dnsmasq</strong> service is very useful in that situation.</p>

[![dnsmasq schema](/assets/img/posts/rpiqa/dnsmasq.png)](/assets/img/posts/rpiqa/dnsmasq.png)
*dnsmasq schema*

<p>The clients connected to the WiFi will resolve the DNS queries thanks to the dnsmasq client listening to incoming connections. This service is able to perform custom DNS resolutions based on a file that works like a /etc/hosts file.</p>

```
192.168.56.1   ubuntu.tecmint.lan
192.168.56.10  centos.tecmint.lan
```

<p id="httpssniffer">Regarding the HTTPS Sniffer, the implemented solution implies having the SSL certificates from the company installed in a reverse proxy that terminates the SSL connection and forwards the HTTP traffic to an internal server that stores the decrypted requests on a cache.</p>

[![https proxy setup](/assets/img/posts/rpiqa/httpsproxy.png)](/assets/img/posts/rpiqa/httpsproxy.png)
*dnsmasq schema*

<p>This cache can be queried from the GUI to inspect the requests. There's an additional implementation that allows to run pcap capture software directly on the network interface to inspect the packets from the UI.</p>
<p>This environment requires a little setup, that is the configuration of the proxy on the computers of the QA team.</p>
<p>Finally, the QA team were able to save a lot of time setting up their scenarios. With this app it's only a matter of clicking buttons instead of executing weird scripts</p>
<h2>The code ...</h2>
<p>Don't blame too much on the code quality: the whole project was implemented few years ago and as a side project on one or two weekends</p>
<p><a href="https://github.com/adriangalera/rpitester">https://github.com/adriangalera/rpitester</a></p>
<p>Here's a video of me presenting that project ;)</p>
<p><iframe src="https://www.youtube.com/embed/aUaEF87pdms?rel=0&amp;showinfo=0" width="560" height="315" frameborder="0" allowfullscreen="allowfullscreen"></iframe></p>
