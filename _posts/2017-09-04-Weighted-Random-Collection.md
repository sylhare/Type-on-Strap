---
layout: post
title: Weighted Random Collection
tags: [collection, random]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview
Let's see how to make weighted random collection.
<!--more-->

# What is Weighted Random Collection?
The weighted random collection is the collection to get the value randomly by weighted value.
Let's say there are some values like `A`, `B`. You want to get these value ramdomly but by applying the specified percentage.

`A:30%`, `B:70%`

When you call `get()`, then this collection will return `A` as 70% and `B` 30% occurrences.

# Implementation
You need to add the value with double type number which is represented the perecentage.

```java
    public void add(double weight, E value) {
        if (weight <= 0) throw new IllegalArgumentException("The weight value is greater than 0.");
        total += weight;
        naviMap.put(total, value);
    }
```

You can get the value corresponding to the percentage randomly.
```java
    public E get() {
        double value = random.nextDouble() * total;
        Map.Entry<Double, E> ceilingEntry = naviMap.ceilingEntry(value);
        if (ceilingEntry != null) return ceilingEntry.getValue();
        else return null;
    }
```

Since we use a random function, we will get different results each time.
```java
    WeightedRandomCollection<String> collection = new WeightedRandomCollection<>();
    collection.add(5, "This is 50%!!!!!");
    collection.add(3, "This is 30%!!!");
    collection.add(1, "This is 10%!");

    IntStream.range(1,10)
            .forEach(i-> System.out.println(collection.get()));
```

Result is
```text
This is 50!!!!!%
This is 50!!!!!%
This is 50!!!!!%
This is 50!!!!!%
This is 10!%
This is 50!!!!!%
This is 50!!!!!%
This is 30!!!%
This is 50!!!!!%
```

Full sorce code is here!
[Github Source Code](https://github.com/han1448/weighted-random-collection)