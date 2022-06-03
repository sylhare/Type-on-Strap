---
layout: post
title: Combine two Java Optionals
description: How to combine two Optional results in Java 9 using Optional.or compared with Java 8
author-id: "galera"
categories: [java9, optional]
tags: [java9, optional]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/combine-optional/featured-image.jpeg"
thumbnail: "assets/img/posts/combine-optional/featured-image.jpeg"
image: "assets/img/posts/combine-optional/featured-image.jpeg"
---

I was implementing a searching algorithm and I had to search for two things. In this article I describe how I implemented combining two Java optionals.

<p><!--more--></p>

The algorithm had to do something similar to:

1. Search inside an array the first item starting by A
2. If none is found, search the first item starting by B
3. If none is found, throw an Exception

Common code:

```java
public abstract class ItemFinder {

    protected abstract String findItem(List<String> items);

    protected Optional<String> findFirstItemStartingWithA(List<String> items) {
        return findFirstItemStartingWith("A", items);
    }

    protected Optional<String> findFirstItemStartingWithB(List<String> items) {
        return findFirstItemStartingWith("B", items);
    }

    private Optional<String> findFirstItemStartingWith(String letter, List<String> items) {
        return items.stream()
            .filter(item -> item.startsWith(letter))
            .findFirst();
    }
}
```

## Java 8

The only way to achieve this with Optionals before Java 9 was to cascade the Optional call:

```java
     public String findItem(List<String> items) {
        return findFirstItemStartingWithA(items)
            .orElseGet(() -> findFirstItemStartingWithB(items)
                .orElseThrow(() -> new RuntimeException("no item found")));
    }
```

The alternative flow is cascaded in the `orElseGet`, which is a little bit hard to read

## Java > 9

In Java 9, the `Optional.or` method is introduced, which enhaces a lot the code for these scenarios. The code above can be rewritten in a much more understandable fashion:

```java
    public String findItem(List<String> items) {
        return findFirstItemStartingWithA(items)
            .or(() -> findFirstItemStartingWithB(items))
            .orElseThrow(() -> new RuntimeException("no item found"));
    }
```