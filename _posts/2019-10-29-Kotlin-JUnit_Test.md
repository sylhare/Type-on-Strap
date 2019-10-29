---
layout: post
title: Kotlin JUnit Test
tags: [kotlin, junit, spring]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview

### Using mockito
```kotlin
import org.mockito.Mockito
import kotlin.test.assertTrue

class KotlinFunctionTest {

    // Create mock
    private val mockService = Mockito.mock(BussinessService::class.java)

    fun test() {
        Mockito.`when`(mockService.addStock(100)).thenReturn(true)
        assertTrue(mockService(100))
        Mockito.verify(mockService).addStock(100)
    }
}
```

### More kotlin style
```xml
        <dependency>
            <groupId>com.nhaarman</groupId>
            <artifactId>mockito-kotlin</artifactId>
            <version>1.6.0</version>
            <scope>test</scope>
        </dependency>
```

```kotlin
import org.mockito.Mockito
import kotlin.test.assertTrue

class KotlinFunctionTest {

    // Create mock
    private val mockService = mock()

    fun test() {
        whenever(mockService.addStock(100)).thenReturn(true)
        assertTrue(mockService(100))
    }
}

```

### Reference
* https://site.mockito.org/