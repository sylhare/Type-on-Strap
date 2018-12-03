---
layout: post
title: "자바스크립트로 OOP 구현하기"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "Prototype 기반 언어인 자바스크립트 : 2.OOP 구현하기"
feature-img: "md/img/thumbnail/java-script-logo.png"    
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap: 
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# Prototype 기반 언어인 자바스크립트 : 2.OOP 구현하기

---

프로토타입 기반 언어인 자바스크립트로 OOP를 사용하는 기법

일반 oop 예를 먼저 적절하게 구성하고 js로 샘플 구성하면 되겠다

2단락에서
  2.1 Function 구문을 통해서 객체를 만드는 법
  2.2 인스턴스를 만드는 방법
특히 var a = {} 처럼 오브젝트형 초기화 하는 개념이 명확해야지

함수를 다루는 개념으로 bind 함수를 꼭 함께 설명하면 좋을 것 같아.


[객체지향 프로그래밍_자바스크립트](https://m.blog.naver.com/PostView.nhn?blogId=love_junim&logNo=220584421589&proxyReferer=https%3A%2F%2Fwww.google.co.kr%2F)

### 3 추상화(Abstract)

#### 3.1 prototype 복제 기법

#### 디자인 패턴(추상 팩토리, abstract factory)

[Javascript Abstract Method with ES6](https://medium.com/@yuribett/javascript-abstract-method-with-es6-5dbea4b00027)

[Abstract Factory](https://www.dofactory.com/javascript/abstract-factory-design-pattern)

[Abstract Classes In Javascript](https://ilikekillnerds.com/2015/06/abstract-classes-in-javascript/)

[자바스크립트에서의 추상화란](http://webclub.tistory.com/137)

[자바스크립트 추상 클래스 구현](http://mohwaproject.tistory.com/entry/%EC%9E%90%EB%B0%94%EC%8A%A4%ED%81%AC%EB%A6%BD%ED%8A%B8-%EC%B6%94%EC%83%81-%ED%81%B4%EB%9E%98%EC%8A%A4-%EA%B5%AC%ED%98%84)

[제로 초 - 추상화 팩토리](https://www.zerocho.com/category/JavaScript/post/57b9692ae492d01700b0b75a)

### 4. 캡슐화(Encapsulation)

캡슐화는 낮은 결합도를 유지할 수 있도록 해주는 객체지향 설계 원리다.

 정보은닉을 통해 높은 응집도와 낮은 결합도를 갖도록 한다.
 
 [Class Member Encapsulation in JavaScript: Data Hiding](https://www.htmlgoodies.com/beyond/javascript/class-member-encapsulation-in-javascript-data-hiding.html)
 
#### 3.1 생성자 함수

[자바스크립트 생성자 함수](http://improver.tistory.com/576)

#### 3.2 정보 은닉  - 낮은 결합,  비공개 프로퍼티, 메소드 비공개 함수 선언

[비공개 프로퍼티 메소드1](http://webclub.tistory.com/312)

[비공개 프로퍼티 메소드 ](http://realmojo.tistory.com/74)

[자바스크립트 강좌 캡슐화 정보은닉 ](https://codingcoding.tistory.com/743)

#### 3.3 클로저 패턴 - 디자인 패턴

[캡슐화 - 클로저](http://webclub.tistory.com/387)

[캡슐화 와 클로져](https://meetup.toast.com/posts/90)


### 5. 일반화 관계(Generalization) - 또 다른 캡슐화

 일반화는 여러 개체들이 가진 공통된 특성을 부각시켜 하나의 개념이나 법칙으로 성립시키는 과정이다.

 일반화 관계는 객체지향 프로그래밍 관점에서는 상속 관계 라 한다.
 따라서 속성이나 기능의 재사용만 강조해서 사용하는 경우가 많다.
 하지만 이는 일반화 관계를 극히 한정되게 바라보는 시각이다.
[자바 oop](https://gmlwjd9405.github.io/2018/07/05/oop-features.html)
 
### 5.1 함수 상속(Inheritance)
### 5.1.1 함수 재정의

### 6 메모리 관리 - 디자인 패턴 플라이급, flyweight

[디자인 패턴(플라이급, flyweight)](https://www.zerocho.com/category/JavaScript/post/57bbb0a3f6f59c170071d2e2)



### 6. 다형성(Polymorphism)

### 6.1 함수 오버라이딩 오버로딩

[제로 초 오버라이딩 오버로딩](https://www.zerocho.com/category/JavaScript/post/59c17a58f40d2800197c65d6)


``` javascript
$.fn.extend({
                bizOption : function($$limit, $$errorMesage, $$editable){

                    var $$selector = this.selector;
                    $$selector = $$selector + ">iframe";

                    return $($$selector).get(0).contentWindow.reSetBizOption($$limit, $$errorMesage, $$editable);
                }

                ,bizInfo : function($$type){

                    var $$selector = this.selector;
                    $$selector = $$selector + ">iframe";

                    return $($$selector).get(0).contentWindow.getBizInfo($$type);
                }
            });
$("#file1").bizOption(1, "사원증 첨부파일은 최대 1개까지 등록 가능합니다.", false);
```