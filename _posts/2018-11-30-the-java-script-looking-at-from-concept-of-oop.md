---
layout: post
title: "나만 몰랐던 자바스크립트의 OOP"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "OOP의 개념으로 바라본 자바스크립트"
feature-img: "md/img/thumbnail/java-script-logo.png"
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap:
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# OOP의 개념으로 바라본 자바스크립트

---

 본 포스팅에선 프로토타입 기반 언어인 자바스크립트와 클래스 기반언어의 차이점과 프로토타입 프로그래밍에서의 객체에 대해 중심적으로 다룰 예정이다. 이번 포스팅엔 크게 세 가지의 학습 목표가 있다.

### 학습 목표

- 프로토타입 기반과 클래스 기반의 객체의 차이
- 프로토타입 프로그래밍의 의미
- 자바스크립트의 객체

---

### 프로토 타입 기반 프로그래밍이란


 여기서 우리는 프로토타입 기반 언어가 OOL를 시초로 두고 있다는걸 알 수 있고 프로토타입 기반 언어인 자바스크립트 또한 OOP의 능력을 지니고 있다는걸 알 수 있다.


프로토타입 기반 프로그래밍은 클래스가 존재하지 않는 객체지향 프로그래밍의 한가지 스타일로, 동작 재사용(behavior reuse, 클래스기반 언어에서는 상속이라고함)은 프로토타입으로서 존재하는 객체를 데코레이팅하는 과정을 통해 수행된다.
프로토타입 기반 언어의 원형적인 예는 David Ungar과 Randall Smith가 개발한 'Self'라는 프로그래밍 언어이다. 그러나 클래스가 없는 프로그래밍 스타일이 최근 인기를 얻으며 성장하였고, 자바스크립트, Cecil, NewtonScript, Io, MOO, REBOL, Kevo, Squeak 등의 언어에서 채택되어 왔다.2

---

프로토타입 기반 프로그래밍은 객체지향 프로그래밍의 한 형태의 갈래로 클래스가 없고, 클래스 기반 언어에서 상속을 사용하는 것과는 다르게, 객체를 원형(프로토타입)으로 하여 복제의 과정을 통하여 객체의 동작 방식을 다시 사용할 수 있다. 프로토타입기반 프로그래밍은 클래스리스(class-less), 프로토타입 지향(prototype-oriented) 혹은 인스턴스 기반(instance-based) 프로그래밍이라고도 한다.

프로토타입 기반 언어의 가장 원조격인 프로그래밍 언어인 셀프는 데이비드 엉거와 랜덜 스미스가 개발했다. 클래스리스 프로그래밍은 최근에 와서 많이 유명해졌는데, 자바스크립트와 모픽 프레임워크를 사용하는 스퀵에 적용되었고, 그 외에 세실, 뉴튼스크립트, 아이오, 무, 리볼, 케보 등에 적용되었다.(#프로토타입 기반 언어 목록 참고.)





 [클래스 기반 & 프로토타입 기반 객체 지향 언어](http://skyul.tistory.com/84)

[전체적인 틀 잡기 ](https://developer.mozilla.org/ko/docs/Web/JavaScript/Introduction_to_Object-Oriented_JavaScript)


### 프로토 타입 기반의 객체



3. 프로토 타입을 가지고 있는 자료형 특징
 3.1 가장 기초적인 Object 부터 int float string 등등
 3.2 확장 프로토 타입 DateTime

[Javascript 기초 - Object prototype 이해하기](http://insanehong.kr/post/javascript-prototype/)

[제로초 - 자바스크립트 변수(Variable), 자료형](https://www.zerocho.com/category/JavaScript/post/57271d6e5aec14515b949b4b)

[제로 초 - 실행 컨테스트](https://www.zerocho.com/category/JavaScript/post/5741d96d094da4986bc950a0)

[제로 초 - 함수 선언](https://www.zerocho.com/category/JavaScript/post/572dcbbd2115c895b0f248fd)

### 클래스 기반과의 비교

[클래스 기반(일반 언어)과 프로토타입(자바스크립트) 기반 비교](http://webclub.tistory.com/162)


### 자바스크립트와  프로토 타입 관계 -  prototype과 `__proto__`, constructor

### 객체 복제와 객체 재사용 - bind ,call apply 3가지 개념


[프로토 타입](https://poiemaweb.com/js-prototype)

`__proto__`  표준이 아니라 prototype 을 쓰는게 합리적

[프로토타입 이해](http://www.nextree.co.kr/p7323/)

[제로 초](https://www.zerocho.com/category/JavaScript/post/573c2acf91575c17008ad2fc)

### this의 차이

[JavaScript에서 '프로토 타입'과 'this'의 사용?](https://code.i-harness.com/ko-kr/q/4be560)

[자바스크립트 this 바인딩 우선순위](http://blog.jeonghwan.net/2017/10/22/js-context-binding.html)

[this의 이해](http://webframeworks.kr/tutorials/translate/explanation-of-this-in-javascript-1/)


프로토타입 기반 언어에 필요한 깊은 키워드들이 들어가야한다고 말한거지

모든 OOP 개념을 설명하면서도 절대로 놓치지 않고 유지해야하는 개념은 자바스크립트는 Prototype based Language라는 점


this 개념을 잘 이해해야해

예제를 찾아보면 클래스 this 객체와 내부에서 선언된 that 객체를 예로 들어서 설명이 있을 텐데  기본적으로는 항상 클래스(객체)로 선언된 내부의 자료를 가르키는 개념인건 동일해

자바스크립트의 전역성 때문에 사용되는 위치에 따라서 다르게 이해될 수 있는 측면이 1차적인 문제이고

자바에서도 this는 현재 인스턴스의 Metadata 가 되는 클래스를 직접 가르키고 있는거지

최상단 인스턴스?? 이런 개념은 아니고

### 마무리

 클래스 기반 언어인 자바로 개발을 해왔던 나는 자바스크립트를 개발 자유도가 높은 언어 또는 함수로 이뤄진 언어로 단정을 짓었다. 또한,  객체는 곧 클래스라고 성급한 일반화의 오류를 범해 클래스 개념이 없는 자바스크립트에선 OOP가 불가능하다고 생각했고 그렇게 알고 사용했다.

때문에, 중복코드로 난발되어진 코드를 보며  '추상화, 상속 등 OOP의 장점을 자바스크립트에서 사용할 순 없을까?'라는 생각이 들었다. 이 생각을 발판삼아 고찰한 결과를 본 포스팅에 작성했다.

본 포스팅을 마무리 지었을땐 '프로토타입 기반 언어인 자바스크립트도 강력한 OOP 능력을 지니고 있다.'는 결론을 다다르게 되었고 성급한 일반화의 오류를 범했던 내 무지에 대해 반성하는 시간이 되었다.
