---
layout: post
title: Singleton Pattern
tags: [java, singleton]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview
The singleton pattern is a software design pattern that restricts the instantiation of a class to one single instance. 
This is useful when exactly one object is needed to coordinate actions across the system.
<!--more-->


### Singleton
```java
    private static final Singleton mySingleton = new Singleton();
    
    private Singleton() {
    }

    public static Singleton getInstance() {
        return mySingleton;
    }
```

The above implementation is so simple.
But `Singleton` is created even though client application might not use it.
It is waste of memory when MySingleton has large amount of resources.

Let's try solve this issue.

### Lazy Initialization
```java
    private static LazyMySingleton mySingleton = null;

    private LazyMySingleton() {
    }

    public static LazyMySingleton getInstance() {
        if (mySingleton == null) {
            mySingleton = new LazyMySingleton();
        }
        return mySingleton;
    }
```
This code looks good. The instance will be created when client calls `getInstance()` method, 
But this code can cause issues in the multi-threaded environment.
The multiple threads can be inside the `if` statement at the same time.
In this time, it is possible to create instant more than one.

Let's find other ways.

### Thread Safe Singleton
```java
    private static ThreadSafeMySingleton mySingleton = null;

    private ThreadSafeMySingleton() {
    }

    public static synchronized ThreadSafeMySingleton getInstance() {
        if (mySingleton == null) mySingleton = new ThreadSafeMySingleton();
        return mySingleton;
    }
```

Above implementation works fine and provides thread-safety but it can reduces the performance when we use `synchronized` keyword. 
Because only one thread can access `getInstance`.

### Double Check Thread Safe Singleton
```java
    private static DoubleCheckThreadSafeMySingleton mySingleton = null;

    private DoubleCheckThreadSafeMySingleton() {
    }

    public static DoubleCheckThreadSafeMySingleton getInstance() {
        if (mySingleton == null) {
            synchronized (DoubleCheckThreadSafeMySingleton.class) {
                if (mySingleton == null) 
                    mySingleton = new DoubleCheckThreadSafeMySingleton();
            }
        }
        return mySingleton;
    }
```
We can double check the `null` value and locate `synchronized` keyword in the second check logic.
Now all thread can access `getInstance()` at the same time.
It's better than ThreadSafeMySingleton in terms of performance.

### Helper Singleton
THere is an another approach to create the Singleton class using an inner static helper class.
```java
    private InnerStaticMySingleton() {
    }

    public static InnerStaticMySingleton getInstance() {
        return SingletonHelper.INSTANCE;
    }

    private static class SingletonHelper {
        private static final InnerStaticMySingleton INSTANCE = new InnerStaticMySingleton();
    }
```

Happy coding :)