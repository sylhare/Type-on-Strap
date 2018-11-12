---
layout: post
title: "OOP의 개념으로 바라본 Prototype 기반 언어인 자바스크립트"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "prototype과 `__proto__`, constructor"
feature-img: "md/img/thumbnail/java-script-oop.png"    
thumbnail: "md/img/thumbnail/java-script-oop.png"
excerpt_separator: <!--more-->
sitemap:
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# 나만 몰랐던 자바 스크립트의 OOP

---

# 나만 몰랐던 자바 스크립트 OOP 지향 프로그래밍


[전체적인 틀 잡기 ](https://developer.mozilla.org/ko/docs/Web/JavaScript/Introduction_to_Object-Oriented_JavaScript)

## [#1] 

[[Java] OOP(객체지향 프로그래밍)의 특징](https://gmlwjd9405.github.io/2018/07/05/oop-features.html)

프로토타입 기반 언어에 필요한 깊은 키워드들이 들어가야한다고 말한거지

모든 OOP 개념을 설명하면서도 절대로 놓치지 않고 유지해야하는 개념은 자바스크립트는 Prototype based Language라는 점

### 1단락 일반 OOP와 다른 점은 프로토타입 기반 언어의 특징
 1-2 단락 프로토타입을 가지고 있는 자료형이 순차적으로 설명
    1-2-1 가장 기초적인 Object 부터 int float string 등등
    1-2-2 확장 프로토 타입 DateTime

###  prototype과 `__proto__`, constructor

`__proto__`  표준이 아니라 prototype 을 쓰는게 합리적

[프로토타입 이해](http://www.nextree.co.kr/p7323/)

[제로 초]((https://www.zerocho.com/category/JavaScript/post/573c2acf91575c17008ad2fc)
)

### this

[this의 이해](http://webframeworks.kr/tutorials/translate/explanation-of-this-in-javascript-1/)

this 개념을 잘 이해해야해

예제를 찾아보면 클래스 this 객체와 내부에서 선언된 that 객체를 예로 들어서 설명이 있을 텐데  기본적으로는 항상 클래스(객체)로 선언된 내부의 자료를 가르키는 개념인건 동일해

자바스크립트의 전역성 때문에 사용되는 위치에 따라서 다르게 이해될 수 있는 측면이 1차적인 문제이고

자바에서도 this는 현재 인스턴스의 Metadata 가 되는 클래스를 직접 가르키고 있는거지

최상단 인스턴스?? 이런 개념은 아니고

### bind call apply 3가지 개념

## [#2] 프로토타입 기반 언어인 자바스크립트로 OOP를 사용하는 기법

일반 oop 예를 먼저 적절하게 구성하고 js로 샘플 구성하면 되겠다

2단락에서
  2.1 Function 구문을 통해서 객체를 만드는 법
  2.2 인스턴스를 만드는 방법
특히 var a = {} 처럼 오브젝트형 초기화 하는 개념이 명확해야지

함수를 다루는 개념으로 bind 함수를 꼭 함께 설명하면 좋을 것 같아.

### 추상화

### 캡슐화
캡슐화는 낮은 결합도를 유지할 수 있도록 해주는 객체지향 설계 원리다.

 정보은닉을 통해 높은 응집도와 낮은 결합도를 갖도록 한다.

 ### 일반화 관계(Generalization)
 일반화는 여러 개체들이 가진 공통된 특성을 부각시켜 하나의 개념이나 법칙으로 성립시키는 과정이다.

 일반화 관계는 객체지향 프로그래밍 관점에서는 상속 관계 라 한다.
 따라서 속성이나 기능의 재사용만 강조해서 사용하는 경우가 많다.
 하지만 이는 일반화 관계를 극히 한정되게 바라보는 시각이다.
 https://gmlwjd9405.github.io/2018/07/05/oop-features.html

 상속


### 함수 상속





### 생성자 함수

### 메소드 정의



### 다형성

### 함수 재정의

### 함수 오버라이딩 오버로딩

[제로 초 오버라이딩 오버로딩](https://www.zerocho.com/category/JavaScript/post/59c17a58f40d2800197c65d6)

### 디자인 패턴 플라이급, flyweight

[디자인 패턴(플라이급, flyweight)](https://www.zerocho.com/category/JavaScript/post/57bbb0a3f6f59c170071d2e2)

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
