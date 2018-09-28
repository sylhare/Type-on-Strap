---
layout: post
title: "TDD : Test Double"
tags: [TDD, TestDouble, Mock, Stub, Fake]
display: "false"
feature-img: "/nyppagesix.files.wordpress.com/2018/04/bobby-hanton-chris-hemsworth.jpg?quality=90&strip=all"              
thumbnail: "/nyppagesix.files.wordpress.com/2018/04/bobby-hanton-chris-hemsworth.jpg?quality=90&strip=all"
subtitle: "위험에 대처하는 자세"
excerpt_separator: <!--more-->
sitemap:
changefreq: daily
priority: 1.0
---



<!--more-->

# 위험에 대처하는 자세

---

### Test double

테스트 더블는 영화에서 나오는 스턴트 대역(Stunt double)에서 유례된 말이다. 테스트시에 실제 객체를 대신 할 수 있는 객체를 의미 합니다. (여기서 ‘더블’이란 말의 유래는 영화에서 스턴트 대역[stunt double]을 생각하시면 될 듯 하네요.)


[https://en.wikipedia.org/wiki/Test_double](https://en.wikipedia.org/wiki/Test_double)

[https://stackoverflow.com/questions/12827580/mocking-vs-spying-in-mocking-frameworks](https://stackoverflow.com/questions/12827580/mocking-vs-spying-in-mocking-frameworks)

[https://eminentstar.github.io/2017/07/24/about-mock-test.html](https://eminentstar.github.io/2017/07/24/about-mock-test.html)
[http://www.jpstory.net/2013/07/26/know-your-test-doubles/](http://www.jpstory.net/2013/07/26/know-your-test-doubles/)

테스트 더블이란 테스트시에 실제 객체를 대신 할 수 있는 객체를 의미 합니다. (여기서 ‘더블’이란 말의 유래는 영화에서 스턴트 대역[stunt double]을 생각하시면 될 듯 하네요.)

테스트 더블을 Mock Object 로 흔히들 알고 있는데요. 그 종류로는 Stub, Mock, Fake Object 등이 있습니다. 각각 다른 용도를 가지고 있기 때문에 어떤 테스트 더블을 언제 써야할지 알기 위해서 서로 구분하는 것은 중요합니다.

각 종류에 대해 간단히 살펴보면,

Stub은 로직이 없고 단지 원하는 값을 반환합니다. 테스트시에 “이 객체는 무조건 이 값을 반환한다”고 가정할 경우 사용할 수 있습니다. Stub은 보통 작성하기 쉽지만 불필요한 boilerplate 코드를 줄이기 위해서 Mocking Framework을 이용하는게 편합니다.

> [Types of test doubles](https://en.wikipedia.org/wiki/Test_double) <br/>
> Test stub <br/>
> Mock object <br/>
> Test spy <br/>
> Fake object <br/>
> Dummy object <br/>

---

### Stub

_used for providing the tested code with "indirect input"_

---

### Mock

_used for verifying "indirect output" of the tested code, by first defining the expectations before the tested code is executed_

---

### Spy

_used for verifying "indirect output" of the tested code, by asserting the expectations afterwards, without having defined the expectations before the tested code is executed. It helps in recording information about the indirect object created_

---

### Fake

_used as a simpler implementation, e.g. using an in-memory database in the tests instead of doing real database access_

---

### Dummy

_used when a parameter is needed for the tested method but without actually needing to use the parameter_

---

### 마무리

---

### 참고

https://adamcod.es/2014/05/15/test-doubles-mock-vs-stub.html

https://blog.pragmatists.com/test-doubles-fakes-mocks-and-stubs-1a7491dfa3da

https://martinfowler.com/bliki/TestDouble.html

https://github.com/testdouble/contributing-tests/wiki/Test-Double

https://lostechies.com/derekgreer/2011/05/15/effective-tests-test-doubles/

https://laurentkempe.com/2010/07/17/Unit-Test-using-test-doubles-aka-Mock-Stub-Fake-Dummy/