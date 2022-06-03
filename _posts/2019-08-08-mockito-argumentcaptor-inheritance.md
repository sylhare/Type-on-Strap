---
layout: post
title: Mockito ArgumentCaptor with inheritance
description: Working with Mockito's ArgumentCaptor I discover it does not work as expected with child classes. This article describe a workaround to keep using it.
author-id: "galera"
categories: [testing]
tags: [testing,java,mockito]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/mockito-argcaptor-inheritance/featured.jpg"
thumbnail: "assets/img/posts/mockito-argcaptor-inheritance/featured.jpg"
image: "assets/img/posts/mockito-argcaptor-inheritance/featured.jpg"
---

Working with Mockito's ArgumentCaptor I discover it has a awful issue with inheritance. 
<p><!--more--></p>
Let's suppose we have a parent class named `Animal` and two child classes `Dog` and `Cat`:

```java
public class Animal {

    private String species;

    public Animal(String species) {
        this.species = species;
    }

    public String getSpecies() {
        return species;
    }
}
public class Cat extends Animal {

    private final String colorEyes;

    public Cat(String colorEyes) {
        super("Cat");
        this.colorEyes = colorEyes;
    }
}
public class Dog extends Animal {

    private String name;

    public Dog(String name) {
        super("Dog");
        this.name = name;
    }

    public String getName() {
        return name;
    }
}
```

And a class that accept an argument of type `Animal`:
```java
public class AnimalProcessor {

    public void processAnimal(Animal animal) {
        System.out.println(animal.getSpecies());
    }
}
```

One might think on writing the unit test of `AnimalProcessor` similar to this snippet:

```java
public class ArgumentCaptorInheritanceTest {

    @Mock
    private AnimalProcessor animalProcessor;

    @Before
    public void init() {
        MockitoAnnotations.initMocks(this);
    }

    @Test
    public void shouldProcessDog() {
        Dog dog = new Dog("Rex");
        Cat cat = new Cat("blue");

        ArgumentCaptor<Dog> dogArgumentCaptor = ArgumentCaptor.forClass(Dog.class);

        animalProcessor.processAnimal(dog);
        animalProcessor.processAnimal(cat);

        Mockito.verify(animalProcessor).processAnimal(dogArgumentCaptor.capture());

        Assert.assertEquals("Rex", dogArgumentCaptor.getValue().getName());
    }
```

Well, this fails ... `ArgumentCaptor` does not work well in this case. It makes the `verify` fail because the method has been called twice. 

The expected behaviour is that only `verify` analyses the calls when the `Dog` instance is passed.

In order to execute the test in this way, some ugly workaround needs to be done:

```java
    @Test
    public void shouldProcessDog() {
        Dog dog = new Dog("Rex");
        Cat cat = new Cat("blue");

        ArgumentCaptor<Animal> animalCaptor = ArgumentCaptor.forClass(Animal.class);

        animalProcessor.processAnimal(dog);
        animalProcessor.processAnimal(cat);

        Mockito.verify(animalProcessor, Mockito.times(2)).processAnimal(animalCaptor.capture());

        List<Animal> processedAnimals = animalCaptor.getAllValues();
        Optional<Animal> dogOptional = processedAnimals.stream()
                                        .filter(a -> a instanceof Dog)
                                        .findFirst();
        Assert.assertTrue(dogOptional.isPresent());
        Assert.assertEquals("Rex", ((Dog) dogOptional.get()).getName());
    }
```

Instead of capturing the arguments for `Dog`, you can do it for `Animal`. 

Then the `verify` will successfully capture the two calls and then all the captured values can be analysed. You can filter the captured values for the objects that are of the interested instance.

This example comes from: 

<a href="https://github.com/adriangalera/java-sandbox/tree/master/src/test/java/mockito/argcaptor">https://github.com/adriangalera/java-sandbox/tree/master/src/test/java/mockito/argcaptor</a>

This is a known issue (already reported as an issue in their repo): 

<a href="https://github.com/mockito/mockito/issues/565">https://github.com/mockito/mockito/issues/565</a>

As of the day of writing the article, the issues is still there and it has been opened from 2016 ...