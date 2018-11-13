---
layout: post
title: "프로토타입 기반 언어인 자바스크립트로 OOP를 사용하는 기법"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "프로토타입 기반 자바스크립트로 OOP 기법"
feature-img: "md/img/thumbnail/java-script-logo.png"    
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap: 
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# 프로토타입 기반 언어인 자바스크립트로 OOP를 사용하는 기법

---
프로토타입 기반 언어인 자바스크립트로 OOP를 사용하는 기법

일반 oop 예를 먼저 적절하게 구성하고 js로 샘플 구성하면 되겠다

2단락에서
  2.1 Function 구문을 통해서 객체를 만드는 법
  2.2 인스턴스를 만드는 방법
특히 var a = {} 처럼 오브젝트형 초기화 하는 개념이 명확해야지

함수를 다루는 개념으로 bind 함수를 꼭 함께 설명하면 좋을 것 같아.

### 1. 생성자 함수

### 2. 메소드 정의

### 3 추상화
#### 3.1 prototype 복제 기법

### 4. 캡슐화
#### 4.1 정보 은닉  - 낮은 결합,  비공개 함수 선언

캡슐화는 낮은 결합도를 유지할 수 있도록 해주는 객체지향 설계 원리다.

 정보은닉을 통해 높은 응집도와 낮은 결합도를 갖도록 한다.

 ### 5. 일반화 관계(Generalization) - 또 다른 캡슐화
 ### 5.1 함수 상속
 ### 5.2 함수 재정의

 일반화는 여러 개체들이 가진 공통된 특성을 부각시켜 하나의 개념이나 법칙으로 성립시키는 과정이다.

 일반화 관계는 객체지향 프로그래밍 관점에서는 상속 관계 라 한다.
 따라서 속성이나 기능의 재사용만 강조해서 사용하는 경우가 많다.
 하지만 이는 일반화 관계를 극히 한정되게 바라보는 시각이다.
 https://gmlwjd9405.github.io/2018/07/05/oop-features.html

### 6. 다형성
### 6.1 함수 오버라이딩 오버로딩

[제로 초 오버라이딩 오버로딩](https://www.zerocho.com/category/JavaScript/post/59c17a58f40d2800197c65d6)

### 7. 메모리 관리 - 디자인 패턴 플라이급, flyweight

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

