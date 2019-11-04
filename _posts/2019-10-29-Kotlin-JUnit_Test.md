---
layout: post
title: Kotlin JUnit Test
tags: [kotlin, junit, spring]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview
Kotlin JUnit5 Test 
### Dependencies
```xml
    <properties>
        <java.version>1.8</java.version>
        <kotlin.version>1.3.50</kotlin.version>
        <mockito-kotlin.version>1.6.0</mockito-kotlin.version>
        <mockito-inline.version>3.1.0</mockito-inline.version>
    </properties>

   <dependencies>
        <dependency>
            <groupId>org.jetbrains.kotlin</groupId>
            <artifactId>kotlin-stdlib-jdk8</artifactId>
            <version>${kotlin.version}</version>
        </dependency>
        <dependency>
            <groupId>org.jetbrains.kotlin</groupId>
            <artifactId>kotlin-test-junit5</artifactId>
            <version>${kotlin.version}</version>
            <scope>test</scope>
        </dependency>
        <dependency>
            <groupId>org.jetbrains.kotlin</groupId>
            <artifactId>kotlin-test-junit</artifactId>
            <version>${kotlin.version}</version>
            <scope>test</scope>
        </dependency>
        <dependency>
            <groupId>com.nhaarman</groupId>
            <artifactId>mockito-kotlin</artifactId>
            <version>${mockito-kotlin.version}</version>
            <scope>test</scope>
        </dependency>
        <dependency>
            <groupId>org.mockito</groupId>
            <artifactId>mockito-inline</artifactId>
            <version>${mockito-inline.version}</version>
            <scope>test</scope>
        </dependency>
    </dependencies>

    <build>
        ...
    </build>
```

### Test class
```kotlin
class CalcService(private val sumService: SumService) {
    fun addSum(a: Int, b: Int): Int = sumService.add(a, b)
}

class SumService {
    fun add(a: Int, b: Int): Int = a + b
}

```

### Using mockito
#### Verify
```kotlin
import org.mockito.Mockito
import kotlin.test.assertTrue

class KotlinFunctionTest {

    private val mockService: SumService = mock()

    @Test
    fun addWithVerify() {
        whenever(mockService.add(1, 1)).thenReturn(2)
        val calcService = CalcService(mockService)
        calcService.addSum(1, 1)

        verify(mockService).add(1, 1)
    }
}
```
#### Assertion
```kotlin
import org.mockito.Mockito
import kotlin.test.assertTrue

class KotlinFunctionTest {

    private val mockService: SumService = mock()

    @Test
    fun addWithAssertion() {
        val expected = 2
        whenever(mockService.add(1, 1)).thenReturn(expected)
        val calcService = CalcService(mockService)

        val actual = calcService.addSum(1, 1)
        assertEquals(expected, actual)
    }
}
```

### Reference
* https://site.mockito.org/