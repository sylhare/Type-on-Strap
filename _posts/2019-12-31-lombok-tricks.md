---
layout: post
title: Lombok tricks
description: Compilation of lombok tricks and useful techniques
author-id: "galera"
categories: [java,lombok]
tags: [java,lombok]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/lombok/featured-image.jpg"
thumbnail: "assets/img/posts/lombok/featured-image.jpg"
image: "assets/img/posts/lombok/featured-image.jpg"
---

<p>
Lombok is this very beatiful tool to reduce the burden of writing Java code, but sometimes it could 
be hard to tame. In this article I write down some issues and solutions I found while using lombok.
</p>
<p><!--more--></p>

<h2>Inmutable objects and Jackson</h2>

Let's say we want to have an inmutable object (`@Value`) such as:

```java
@Value
@Builder
public class Foo {

    private String id;
    private String description;
}
```

If that's the structure returned by some API, one could do the following to consume it:

```java
RestTemplate restTemplate = new RestTemplate();
HttpEntity<String> entity = new HttpEntity<>();
try {
    ResponseEntity<Foo> response = restTemplate.exchange(url HttpMethod.GET, entity,Foo.class);
    return Optional.ofNullable(response.getBody());
} catch (Exception ex) {
    log.error("Error requesting to API: {}", ex);
}
```

However, in this case, some exception like the following will be thrown by Jackson JSON deserialization library:

```
Caused by: com.fasterxml.jackson.databind.exc.InvalidDefinitionException: 
Cannot construct instance of `com.gal.Foo` (no Creators, like 
default construct, exist): cannot deserialize from Object value 
(no delegate- or property-based Creator)
```

The solution for this is pretty simple, we need to configure lombok with a private no-arg constructor and a constructor with all arguments:

```java
@Value
@NoArgsConstructor(force = true, access = AccessLevel.PRIVATE)
@AllArgsConstructor
@Builder
public class Foo {

    private String id;
    private String description;
}
```

This way, Jackson can deserialize the object with minimal lombok configuration