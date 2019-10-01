---
layout: post
title: Sprig Redis Json Example
tags: [Spring, Redis]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview
This is simple example of Spring-Boot-Redis.
<!--more-->

# Version
* java8
* Redis : 3.2.8


# Pom.xml
Add spring boot dependencies
```xml
<parent>
    <groupId>org.springframework.boot</groupId>
    <artifactId>spring-boot-starter-parent</artifactId>
    <version>1.5.2.RELEASE</version>
</parent>

<properties>
    <project.build.sourceEncoding>UTF-8</project.build.sourceEncoding>
</properties>

<dependencies>
    <dependency>
        <groupId>org.springframework.boot</groupId>
        <artifactId>spring-boot-starter-data-redis</artifactId>
    </dependency>
    <!-- https://mvnrepository.com/artifact/com.fasterxml.jackson.core/jackson-core -->
    <dependency>
        <groupId>com.fasterxml.jackson.core</groupId>
        <artifactId>jackson-core</artifactId>
        <version>2.8.7</version>
    </dependency>
    <!-- https://mvnrepository.com/artifact/com.fasterxml.jackson.core/jackson-databind -->
    <dependency>
        <groupId>com.fasterxml.jackson.core</groupId>
        <artifactId>jackson-databind</artifactId>
        <version>[2.8.11.1,)</version>
    </dependency>
    <!-- https://mvnrepository.com/artifact/org.apache.commons/commons-lang3 -->
    <dependency>
        <groupId>org.apache.commons</groupId>
        <artifactId>commons-lang3</artifactId>
        <version>3.5</version>
    </dependency>

    <dependency>
        <groupId>junit</groupId>
        <artifactId>junit</artifactId>
        <scope>test</scope>
    </dependency>
</dependencies>
```

# Create AppConfig.java
Create `@Bean` to use `RedisTemplate`.
Add `GenericJackson2JsonRedisSerializer` to convert to json and set valueSerializer.

```java
@Bean
public <T> RedisTemplate<String, T> redisTemplate(RedisConnectionFactory connectionFactory) {
	RedisTemplate<String, T> template = new RedisTemplate<String, T>();
	template.setConnectionFactory(connectionFactory);
	template.setKeySerializer(new StringRedisSerializer());
	template.setValueSerializer(new GenericJackson2JsonRedisSerializer());
	return template;
}
```

# Create model
Create model class to save in redis.
```java
public class Author implements Serializable {
	int id;
	String name;
	int age;
	List<Book> books;
```

```java
public class Book implements Serializable {
	int bookId;
	String title;
	String description;
	Date publishDate;
```

# Create application.yml
Create application.yml which has information of redis connection.
```yaml
spring:
  redis:
    host: localhost
    port: 6379
```


# Create main class
Add template instance by using @Autowired
```java
@Autowired
private RedisTemplate<String, Author> template;
```

Add main function to execute spring boot application
```java
public static void main(String[] args) {
	SpringApplication.run(App.class, args);
}
```

To execute this app, add implements CommandLineRunner.
```java
@SpringBootApplication
public class App implements CommandLineRunner {
    @Override
    public void run(String... arg0) throws Exception {
        
    }
}
```
And override run() method.

This `run()` method will insert author into redis and query again.
```java
System.out.println("start");
Author author = getAuthor();
ValueOperations<String, Author> ops = this.template.opsForValue();
String key = "book:name:" + author.getId();
System.out.println(key);
ops.set(key, author);
Author fromDbAuthor;
fromDbAuthor = ops.get(key);
System.out.println(fromDbAuthor);
```

After execute this app, you can see below logs.
```text
...
2017-03-06 22:20:30.289 INFO 802 --- [ main] o.s.j.e.a.AnnotationMBeanExporter : Registering beans for JMX exposure on startup 
book:name:1 
com.oppalove.spring.redis.model.Author@3dd69f5a[id=1,name=Han SeungHyeon,age=38,
books=[com.oppalove.spring.redis.model.Book@59a67c3a[bookId=1,title=Gogo land,description=i don't know,publishDate=<null>], 
com.oppalove.spring.redis.model.Book@5003041b[bookId=2,title=Java world,description=Let's study,publishDate=<null>]]] 
2017-03-06 22:20:30.568 INFO 802 --- [ main] com.oppalove.spring.redis.App : Started App in 3.664 seconds (JVM running for 4.502)
...
```

# Check Redis
Connect to redis database
```text
seungui-MacBook-Air:redis-3.2.8 seunghyunhan$ ./src/redis-cli 
127.0.0.1:6379>
```

Query by using key book:name:1
```text
127.0.0.1:6379> get book:name:1 
"{"@class":"com.oppalove.spring.redis.model.Author","id":1,"name":"Han SeungHyeon","age":38,
"books":["java.util.ArrayList",[{"@class":"com.oppalove.spring.redis.model.Book","bookId":1,"title":"Gogo land","description":"i don't know","publishDate":null},
{"@class":"com.oppalove.spring.redis.model.Book","bookId":2,"title":"Java world","description":"Let's study","publishDate":null}]]}" 
127.0.0.1:6379>
```

[Github Source Code](https://github.com/han1448/spring-boot-redis-example)