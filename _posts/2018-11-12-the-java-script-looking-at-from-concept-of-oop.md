---
layout: post
title: "OOP의 개념으로 바라본 자바스크립트"
tags: [JavaScript, OOP]
categories: [JavaScript]
subtitle: "프로토 타입 기반 언어인 자바스크립트의 OOP"
feature-img: "md/img/thumbnail/java-script-logo.png"
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap:
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# 프로토타입 기반 언어인 자바스크립트의 OOP

---

 일반적으로 객체 지향 언어는 클래스 기반 언어와 프로토타입 기반 언어로 나눌수 있다.

- 객체 지향 언어
	- 클래스 기반언어 - Java
	- 프로토타입 기반언어 - JavaScript

흔히 OOP의 대표 언어인 JAVA

---

[OOP 객체지향 프로그래밍](http://www.terms.co.kr/OOP.htm)


객체지향 프로그래밍(이하 줄여서 'OOP'라 칭함)은 컴퓨터 프로그램의 개발을 완전히 새로운 시각으로 바라다보는 혁명적 개념이라 할 수 있는데, 동작보다는 객체, 논리보다는 자료를 바탕으로 구성된다. 프로그램은 전통적으로 논리적인 수행 즉, 입력을 받아 처리한 다음, 결과를 내는 것이라는 생각이 지배적이었다. 또한 프로그래밍을 한다는 것은 어떻게 자료를 정의할까 보다는 어떻게 논리를 써나가는 것인가로 간주되었다.

그러나 OOP는 프로그램에서 정말 중요한 것이 논리보다는 오히려 다루고자 하는 객체라는 시각에서 접근하고 있다. 객체의 예로는, 사람(이름, 주소 등으로 묘사되는)에서부터 건물까지, 그리고 상품 판매를 위한 매장(특성이 서술되고 다뤄질 수 있는)에서부터 컴퓨터 바탕화면의 아주 작은 요소인 버튼이나 스크롤바 같은 것들까지를 모두 망라한다.

OOP에서의 첫 단계는 다루고자 하는 모든 객체와, 그것들이 서로 어떤 연관성이 있는지를 식별하는 - 흔히 데이터 모델링이라고 부르는 - 작업이다. 일단 모든 객체를 식별했으면, 객체 클래스로 일반화하고 (플라톤의 "이상적" 의자 개념이 모든 의자를 대표한다고 생각하는 식으로), 그것이 담고 있는 데이터의 종류와 그것을 다룰 수 있는 모든 논리 순서를 정의한다.

논리 순서는 메쏘드라고 부르며, 클래스의 실제 인스턴스를 하나의 "객체"라 하거나, 어떤 상황에서는 하나의 "클래스 활성체"라 한다. 객체 또는 활성체는 컴퓨터 내에서 실제로 수행되는 것이다. 메쏘드는 컴퓨터 명령어를 규정하고, 클래스 객체의 특성은 관련 데이터를 규정한다.
 
OOP에 사용된 개념과 규칙은 다음과 같은 중요한 이점을 제공한다.

데이터 클래스 개념은 일부 또는 모든 주 클래스 특성을 공유하는 데이터 객체의 부 클래스를 정의할 수 있게 한다. 상속이라 불리는 이 OOP 특성은 면밀한 자료 분석과 개발시간 단축, 그리고 좀더 정확한 코딩을 보증하는 효과가 있다.

클래스는 단지 관련된 데이터만 정의하기 때문에, 그 클래스의 인스턴스가 수행될 때 다른 프로그램에서 사용하는 데이터를 (우연이라도) 건드릴 수 없게 된다. 이런 자료 숨김 특성은 높은 시스템 보안을 제공하고, 의도하지 않은 자료의 훼손을 방지한다.

클래스의 정의는 최초로 생성한 프로그램뿐 아니라 다른 OOP에 의해 재사용될 수 있다 (그리고, 이런 이유로 네트웍에 쉽게 분산 사용된다).

데이터 클래스 개념은 언어에 정의되지 않은 새로운 데이터 형식을 프로그래머가 임의로 정의할 수 있게 한다.
Smalltalk는 최초의 객체지향 프로그래밍 언어 중 하나이며, C++와 Java는 최근 가장 인기있는 객체지향 프로그래밍 언어이다. C++의 부분집합이라고 할수 있는 Java는 특히 기업이나 인터넷의 분산 응용프로그램에 사용되도록 설계되었다.

---

  프로토타입 기반 프로그래밍은 OOP에서 파생된 하나의 기법이다.

Prototype 기반 언어인 자바스크립트는 강력한 OOP 능력을 지니고 있다.

본 포스팅에서

[객체 지향 언어의 두 가지 줄기](http://mohwa.github.io/blog/javascript/2015/10/16/prototype/)

[클래스 기반 & 프로토타입 기반 객체 지향 언어](http://skyul.tistory.com/84)


---

 객체지향 프로그래밍은 실제 세계에 기반한 모델을 만들기 위해 추상화를 사용하는 프로그래밍 패러다임이다.


객체지향 프로그래밍은 보다 유연하고 유지보수성이 높은 프로그래밍을 하도록 의도되었고, 대규모 소프트웨어 공학에서 널리 알려져 있다. 객체지향 프로그래밍이 갖는 modularity에 기반한 강력한 힘에 의해, 객체지향적인 코드는 개발을 보다 단순하게 했고, 시간이 흐른 뒤에도 보다 쉽게 이해할 수 있도록 했으며, 복잡한 상황이나 절차들을 덜 모듈화된 프로그래밍 방법들보다 더 직접적으로 분석하고, 코딩하고, 이해할 수 있도록 만들었다.2



[전체적인 틀 잡기 ](https://developer.mozilla.org/ko/docs/Web/JavaScript/Introduction_to_Object-Oriented_JavaScript)

---

### 프로토 타입 기반 프로그래밍이란

프로토타입 기반 프로그래밍은 클래스가 존재하지 않는 객체지향 프로그래밍의 한가지 스타일로, 동작 재사용(behavior reuse, 클래스기반 언어에서는 상속이라고함)은 프로토타입으로서 존재하는 객체를 데코레이팅하는 과정을 통해 수행된다.
프로토타입 기반 언어의 원형적인 예는 David Ungar과 Randall Smith가 개발한 'Self'라는 프로그래밍 언어이다. 그러나 클래스가 없는 프로그래밍 스타일이 최근 인기를 얻으며 성장하였고, 자바스크립트, Cecil, NewtonScript, Io, MOO, REBOL, Kevo, Squeak 등의 언어에서 채택되어 왔다.2

---

프로토타입 기반 프로그래밍은 객체지향 프로그래밍의 한 형태의 갈래로 클래스가 없고, 클래스 기반 언어에서 상속을 사용하는 것과는 다르게, 객체를 원형(프로토타입)으로 하여 복제의 과정을 통하여 객체의 동작 방식을 다시 사용할 수 있다. 프로토타입기반 프로그래밍은 클래스리스(class-less), 프로토타입 지향(prototype-oriented) 혹은 인스턴스 기반(instance-based) 프로그래밍이라고도 한다.

프로토타입 기반 언어의 가장 원조격인 프로그래밍 언어인 셀프는 데이비드 엉거와 랜덜 스미스가 개발했다. 클래스리스 프로그래밍은 최근에 와서 많이 유명해졌는데, 자바스크립트와 모픽 프레임워크를 사용하는 스퀵에 적용되었고, 그 외에 세실, 뉴튼스크립트, 아이오, 무, 리볼, 케보 등에 적용되었다.(#프로토타입 기반 언어 목록 참고.)

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

[[Java] OOP(객체지향 프로그래밍)의 특징](https://gmlwjd9405.github.io/2018/07/05/oop-features.html)

프로토타입 기반 언어에 필요한 깊은 키워드들이 들어가야한다고 말한거지

모든 OOP 개념을 설명하면서도 절대로 놓치지 않고 유지해야하는 개념은 자바스크립트는 Prototype based Language라는 점


this 개념을 잘 이해해야해

예제를 찾아보면 클래스 this 객체와 내부에서 선언된 that 객체를 예로 들어서 설명이 있을 텐데  기본적으로는 항상 클래스(객체)로 선언된 내부의 자료를 가르키는 개념인건 동일해

자바스크립트의 전역성 때문에 사용되는 위치에 따라서 다르게 이해될 수 있는 측면이 1차적인 문제이고

자바에서도 this는 현재 인스턴스의 Metadata 가 되는 클래스를 직접 가르키고 있는거지

최상단 인스턴스?? 이런 개념은 아니고
```
