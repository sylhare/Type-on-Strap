---
layout: post
title: "TDD : Test in Spring"
tags: [TDD]
subtitle: "Chapter3 Spring에서 테스트"
excerpt_separator: <!--more-->
display: "false"
sitemap:
changefreq: daily
priority: 1.0
---

<!--more-->

# Chapter3 Spring에서 테스트

---

### Spring-test를 사용해야 되는 이유

 최근 소프트웨어의 테스트에 대한 책임이 개발자에게 주어지는 추세이고 IT 회사들이 `TDD`의 경험이 있는 개발자들을 선호하고 있다.<br/>
 이에 개발자들이 소프트웨어의 테스트에 대한 관심이 많아지고 있고 프로젝트에서 `TDD`를 도입하려는 행동이 보여지고 있다. 이를 입증하듯 `OKKY` 개발자 커뮤니티에서 `TDD`를 주제로 한 세미나가 최단시간 마감이 되었고 켄트백의 `Test Driven Development`라던가 로버트 C. 마틴의 `Clean Code`가 필독서로 추천하는 이유이다.



많은 개발자들이 스프링 테스트를 사용하고 있고 테스트에 대한 개발자의  
꾸준히 TDD가 각광 받으면서 많은 프로젝트에서 TDD를 도입하기 때문에 사용을 하는걸까 ? <br/>
반적인 개발과정을 보면 Spring-test를 사용해야 되는 이유를 알 수 있다.

#### 일반적인 개발과정

#### Spring-test 개발과정 

---

### 스프링 MVC 테스트 범위

<img src="/md/img/test-in-spring/spring-mvc-architecture.png" height="500px">
<em>Spring MVC 아키텍처</em>

스프링에서의 테스트의 범위는 `Dispatcher Servlet`, `HandlerMapping`, `ViewResolver`, `Controller` 부분이다. 막상 테스트를 진행하려 보니 `Dispatcher Servlet`, `HandlerMapping`, `ViewResolver`은 추상적인 개념이라 테스트를 하기 어려울 수 있다. 하지만 추상적인 개념들은 스프링 테스트에서 처리를 해주기 때문에 개발자는 `Controller`에 해당되는 비즈니스 로직 테스트에 집중하면 된다.

---

### Spring-test 종류


---

### 통합 테스트(Integration Tests)

---

 
 스프링에서의 테스트란 기존 
 
Controller, Service, Model 디자인 패턴으로 설계되어 있어 각 계층의 성격에 맞는 테스트 케이스를 만들어야 된다.