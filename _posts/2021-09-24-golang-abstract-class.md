---
layout: post
title: Golang abstract class
description: How to implement an abstract class on golang
author-id: "galera"
categories: [golang, design-patterns, architecture]
tags: [golang, design-patterns, architecture]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/golang-abstract/featured-image.jpg"
thumbnail: "assets/img/posts/golang-abstract/featured-image.jpg"
image: "assets/img/posts/golang-abstract/featured-image.jpg"
---

I'm doing a side project using golang, and I have a use case where I'd use an abstract class in Java.  Unfortunately, in golang the concept of classes does not exist. 

In this article I describe how can I implement the behaviour I want without abstract class.

<p><!--more--></p>

I'm implementing an alert system in golang. When the alert needs to be activated, I want to play a sound through speakers and blink some LEDs for some period of time.

The behavior is quite simple, it needs to provide an implementation to enable the signaler, disable and query.

We will create two clases, one for playing sounds and another one for blinking the LEDs. And an abstract class to implement the shutdown after some period of time. Something similar to:

```java
abstract class Signaler {
    public void enableForTime(int seconds) {
        //Implementation to enable the signaler for x seconds
    }
    public void disable() {
        //Disable the signalers and cancel any pending timer
    }
    abstract void enable();
    abstract void disableSignaler();
}
class LedSignaler extends Signaler {
    @Override
    void enable() {
        //enable the LEDs
    }
    @Override
    void disableSignaler() {
        //disable the LEDs
    }
}
class SoundSignaler extends Signaler {
    @Override
    void enable() {
        //enable the speaker
    }
    @Override
    void disableSignaler() {
        //disable the speaker
    }
}
```

The idea of this is to abstract the concrete implementation to the caller of the signaler. 

This representation will work in any language that supports inheritance and abstract classes.

## Decorator pattern

In golang, the concept of classes does not exist. So, we need to re-architecture the pattern. 

Instead of using the abstract class, we can re-think the implementation to use a decorator pattern: <a href="https://en.wikipedia.org/wiki/Decorator_pattern">https://en.wikipedia.org/wiki/Decorator_pattern</a>.

> The decorator pattern is a design pattern that allows behavior to be added to an individual object, dynamically, without affecting the behavior of other objects from the same class

It works by defining an interface that will have multiple implementation, we can add new functionalities by adding new implementations. Let's go define the interface of the signaler:

```golang
type Signaler interface {
	Enable()
	IsEnabled() bool
	Disable()
}
```
Now, we can do the implementation of the sound signaler. To implement and interface in golang you need to provide a struct that has the same methods as the interface.

```golang
type soundSignaler struct {
    
}
func (s *soundSignaler) Enable() {
	s.enabled = true
    //start playing sound
}

func (s *soundSignaler) Disable() {
	s.enabled = false
    //stop playing sound
}

func (s *soundSignaler) IsEnabled() bool {
	return s.enabled
}
```
OK, now let's do the implementation of the temporal execution of the signaler. In order to do so, let's create a new type that implement `Signaler`:

```golang
type temporalSignaler struct {
	config          *config.AlarmConfig
	timers          []signalerTimer
	delayedExecutor delayedExecutor
	signaler        Signaler
}

func (s *temporalSignaler) Enable() {
	s.signaler.Enable()
	var t = s.delayedExecutor.executeAfterSeconds(s.config.SecondsStopSignals, func() {
		s.Disable()
	})
	s.timers = append(s.timers, t)
}

func (s *temporalSignaler) Disable() {
	for i := 0; i < len(s.timers); i++ {
		s.timers[i].Stop()
	}
	s.signaler.Disable()
}

func (s *temporalSignaler) IsEnabled() bool {
	return s.signaler.IsEnabled()
}
```
The key point of this struct, is that we're holding a `Signaler` instance on it. We are "decorating" that instance with the temporal disabling functionality. In order to do so, we just need to implement the logic in the type and call the `Signaler` methods of the decorated instance.

The most elegant part of this implementation is that we're defining all the structs as private (name begins with lower case). So, the external modules cannot instantiate them. We can abstract the creation of the sound signaler by creating a factory:

```golang
func NewSoundSignaler(config *config.Config) Signaler {
	return &temporalSignaler{
		delayedExecutor: &defaultDelayedExecutor{},
		config:          &config.Alarm,
		signaler: &soundSignaler{
			config:   config,
			executor: &unixCommandExecutor{},
		},
	}
}
```

External modules can call the factory method, and they will receive a `Signaler` instance, hiding effectively the implementation details of the temporal execution and sound playback.
