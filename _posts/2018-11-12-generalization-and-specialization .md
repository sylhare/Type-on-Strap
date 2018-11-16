---
layout: post
title: "generalization와 specialization "
tags: [OO,  OOL, Class-based language, Prototype-based language]
categories: [OO]
subtitle: "클래스기반 언어와 프로토타입 기반 언어"
feature-img: "md/img/thumbnail/java-script-logo.png"
thumbnail: "md/img/thumbnail/java-script-logo.png"
excerpt_separator: <!--more-->
sitemap:
display: "false"
changefreq: daily
priority: 1.0
---

<!--more-->

# OOL: 클래스기반 언어와 프로토타입 기반 언어

---

 본 포스팅에선 객체지향의 의미와 더불어 OOP에 대해 다룰 예정이다. 이번 포스팅엔 크게 세 가지의 학습 목표가 있다.

### 학습 목표

- OOL의 종류와 개념

###  OOL : Class-based programming, Prototype-based

 객체지향 언어(Object-oriented language, OOL)에는 크게 두가지 언어로 분류된다.

 - OOL
 	- 클래스 기반언어 :  JAVA, C++, C# 등
 	- 프로토타입 기반언어 : JavaScript, ActionScript, JScript 등

  첫 번째가 C++, 자바, C# 등이 공통적으로 쓰고 있는 클래스 기반의 객체지향 언어이고, 또 다른 방식으로 자바 스크립트나 Io 언어 등이 사용하는 프로토타입 방식의 객체 지향 언어가 있습니다.

#### 클래스 기반의 OOL

  클래스 기반 OOL의 경우 클래스가 객체를 찍어내는 틀 역할을 합니다. 객체는 클래스라는 틀에서 붕어빵을 찍듯이 하나씩 찍어내는 개념이죠. C++, 자바 등에 익숙한 사람이라면 new 연산자를 이용해서 클래스로부터 객체를 생성하는 과정을 이해하실 겁니다.

#### 프로토타입 기반의 OOL

  반면에 프로토타입 기반의 언어는 이런 개념을 차용하지 않았습니다. 클래스라는 틀의 존재를 상정치 않고(classless object model이라고도 불립니다), 곧바로 시스템에 여러 종류의 객체가 존재하게 됩니다. 프로토타입 기반의 객체 지향 언어는 이렇게 이미 존재하는 객체를 직접 사용하거나, 클론(clone)해서 사용합니다. 필요에 따라 클론한 객체에 추가적인 기능을 구현(일종의 inheritance)을 하기도 합니다.

  현재는 클래스 기반의 OO가 대세이기 때문에 "객체 지향 = 클래스 기반"이라는 등식이 생겼지만, 실제로는 프로토타입 기반의 OO도 여러 언어에서 사용되고 있습니다. 현재 개발자들에게 가장 익숙한 프로토타입 OO 언어는 아마 JavaScript일 것입니다.

  University of Washington의 Washington Advanced Systems for Programming에 있는 Craig Chambers 교수가 연구했던 프로토타입 기반 객체 지향 언어인 Cecil 프로젝트를 보시면, 프로토타입 기반 OO로 우리가 알고 있는 OO의 기본 개념들을 어떻게 구현하는지 힌트를 얻을 수 있으실 겁니다. 다만 현재 국내 개발자들 사이에서 인지도가 높은 Io 언어는 어떤 차별성이 있는지 잘 모르겠네요.

  반면에 프로토타입 기반의 언어는 이런 개념을 차용하지 않았습니다. 클래스라는 틀의 존재를 상정치 않고(classless object model이라고도 불립니다), 곧바로 시스템에 여러 종류의 객체가 존재하게 됩니다. 프로토타입 기반의 객체 지향 언어는 이렇게 이미 존재하는 객체를 직접 사용하거나, 클론(clone)해서 사용합니다. 필요에 따라 클론한 객체에 추가적인 기능을 구현(일종의 inheritance)을 하기도 합니다.

  현재는 클래스 기반의 OO가 대세이기 때문에 "객체 지향 = 클래스 기반"이라는 등식이 생겼지만, 실제로는 프로토타입 기반의 OO도 여러 언어에서 사용되고 있습니다. 현재 개발자들에게 가장 익숙한 프로토타입 OO 언어는 아마 JavaScript일 것입니다.

  University of Washington의 Washington Advanced Systems for Programming에 있는 Craig Chambers 교수가 연구했던 프로토타입 기반 객체 지향 언어인 Cecil 프로젝트를 보시면, 프로토타입 기반 OO로 우리가 알고 있는 OO의 기본 개념들을 어떻게 구현하는지 힌트를 얻을 수 있으실 겁니다. 다만 현재 국내 개발자들 사이에서 인지도가 높은 Io 언어는 어떤 차별성이 있는지 잘 모르겠네요.

  ---

### 마무리

OOP 혹은 객체지향 프로그래밍이라는 단어는 개발자라면 수도 없이 들었던 단어이다. 특히나 신입 개발자 대상으로  자주 나오는 질문 중 하나로

 클래스 기반 언어인 자바로 개발을 해왔던 나는 자바스크립트를 개발 자유도가 높은 언어 또는 함수로 이뤄진 언어로 단정을 짓었다. 또한,  객체는 곧 클래스라고 성급한 일반화의 오류를 범해 클래스 개념이 없는 자바스크립트에선 OOP가 불가능하다고 생각했고 그렇게 알고 사용했다.

때문에, 중복코드로 난발되어진 코드를 보며  '추상화, 상속 등 OOP의 장점을 자바스크립트에서 사용할 순 없을까?'라는 생각이 들었다. 이 생각을 발판삼아 고찰한 결과를 본 포스팅에 작성했다.

본 포스팅을 마무리 지었을땐 '프로토타입 기반 언어인 자바스크립트도 강력한 OOP 능력을 지니고 있다.'는 결론을 다다르게 되었고 성급한 일반화의 오류를 범했던 내 무지에 대해 반성하는 시간이 되었다.

---

### 참고


 [객체 지향 언어의 두 가지 줄기](http://mohwa.github.io/blog/javascript/2015/10/16/prototype/)
