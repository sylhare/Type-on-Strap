---
layout: post
title: "OOP의 개념으로 바라본 Prototype 기반 언어인 자바스크립트"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "OOP와 프로토타입 기반 언어의 차이점 및 자바스크립트의 동작원리"
feature-img: "md/img/thumbnail/java-script-logo.png"    
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap: 
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# OOP의 개념으로 바라본 Prototype 기반 언어인 자바스크립트

---

[전체적인 틀 잡기 ](https://developer.mozilla.org/ko/docs/Web/JavaScript/Introduction_to_Object-Oriented_JavaScript)

## OOP의 개념으로 바라본 Prototype 기반 언어인 자바스크립트

### 1. 프로토 타입 기반 프로그래밍

### 2. 클래스 기반과의 비교
[클래스 기반(일반 언어)과 프로토타입(자바스크립트) 기반 비교](http://webclub.tistory.com/162)

### 3. 프로토 타입을 가지고 있는 자료형 특징
####   3.1 가장 기초적인 Object 부터 int float string 등등
####   3-2 확장 프로토 타입 DateTime
프로토 타입 복제에 기반을 둘 예정

### 4 객체 복제와 객체 재사용 - bind ,call apply 3가지 개념

### 5. 자바스크립트와  프로토 타입 관계 -  prototype과 `__proto__`, constructor


`__proto__`  표준이 아니라 prototype 을 쓰는게 합리적

[프로토타입 이해](http://www.nextree.co.kr/p7323/)

[제로 초](https://www.zerocho.com/category/JavaScript/post/573c2acf91575c17008ad2fc)

### 6. this의 차이

[JavaScript에서 '프로토 타입'과 'this'의 사용?](https://code.i-harness.com/ko-kr/q/4be560)

[자바스크립트 this 바인딩 우선순위](http://blog.jeonghwan.net/2017/10/22/js-context-binding.html)

[this의 이해](http://webframeworks.kr/tutorials/translate/explanation-of-this-in-javascript-1/)

[[Java] OOP(객체지향 프로그래밍)의 특징](https://gmlwjd9405.github.io/2018/07/05/oop-features.html)

프로토타입 기반 언어에 필요한 깊은 키워드들이 들어가야한다고 말한거지

모든 OOP 개념을 설명하면서도 절대로 놓치지 않고 유지해야하는 개념은 자바스크립트는 Prototype based Language라는 점


this 개념을 잘 이해해야해

예제를 찾아보면 클래스 this 객체와 내부에서 선언된 that 객체를 예로 들어서 설명이 있을 텐데  기본적으로는 항상 클래스(객체)로 선언된 내부의 자료를 가르키는 개념인건 동일해

자바스크립트의 전역성 때문에 사용되는 위치에 따라서 다르게 이해될 수 있는 측면이 1차적인 문제이고

자바에서도 this는 현재 인스턴스의 Metadata 가 되는 클래스를 직접 가르키고 있는거지

최상단 인스턴스?? 이런 개념은 아니고
```
