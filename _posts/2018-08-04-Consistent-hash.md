---
layout: post
title: Consistent Hashing
tags: [Consistent Hash]
author-id: oppalove
excerpt_separator: <!--more-->
---
# What is the Consistent Hashing?
Consistent hashing is a special kind of hashing such that when a hash table is resized, only K/n keys need to be remapped on average, where K is the number of keys, and n is the number of slots. 
In contrast, in most traditional hash tables, a change in the number of array slots causes nearly all keys to be remapped because the mapping between the keys and the slots is defined by a modular operation.

<!--more-->

# Current Hashing Model

The existing hash structure changes the index of existing data if bucket size is increased or decreased.
Below example shows A will store Node B because the hash index points to Node B. B also will store Node C.

![Current Hash Model]({{ "/assets/img/post/2019-09-01/drawit-diagram-14.png" | relative_url}})

But if we add a new Node D, the hash index of the existing nodes will be changed. Whenever a node is added / deleted, all the indexes will be changed relaged to hashing. This will be caused all of data need to be re-partitioned.

![Adding New Node]({{ "/assets/img/post/2019-09-01/drawit-diagram-15.png" | relative_url}})


# Consistent Hashing


Consistent hashing stores data in the closest index node.

![Adding New Node]({{ "/assets/img/post/2019-09-01/drawit-diagram-18.png" | relative_url}})

Data changes on existing nodes are minimized even if new nodes are added. This reduces the re-cache or re-patitioning cost of each server. When a new data is added, it will be stored in the closest node.

![Adding New Node]({{ "/assets/img/post/2019-09-01/drawit-diagram-20.png" | relative_url}})

If NodeD is deleted, data D will be stored in the closest NodeA .

![Adding New Node]({{ "/assets/img/post/2019-09-01/drawit-diagram-22.png" | relative_url}})

The data of each node may not enter uniformly. In this case, a new virtual node is added so as to be distributed as equitably as possible. This allows the nodes to load-balance each node when adding or deleting them.

![Adding New Node]({{ "/assets/img/post/2019-09-01/drawit-diagram-19.png" | relative_url}})

However, it is important to set the maximum number of nodes (including virtual nodes), the number of nodes, and the number of virtual nodes to create evenly distributed nodes.

By using consistent hashing to distribute keys between the servers,  the impact on origin servers will be minimized, preventing potential downtime or performance issues.

[Github Source Code](https://github.com/han1448/consistent_hashing)