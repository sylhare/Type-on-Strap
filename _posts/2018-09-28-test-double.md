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

 프로그램 기능을 테스트시에 해당된 객체를 직접적으로 테스트하기엔 위험요소가 많다. 그 위험요소는 다음과 같다.
 
<img src="/md/img/test-double/test-double_used_reason_of_in_test-code.png">
<em>Test 위험요소</em>

테스트 대상이 될 객체를 직접적으로 사용하면 해당 객체에 사용되는 객체간의 상호관계에서 발생하는 부작용을 야기할 수 있다. 이 때문에 일반적으로 테스트 목적으로 실제 오브젝트와 동작이 같은 주변과 분리 된 독립적인 테스트 객체를 사용한다.

이렇게하면 복잡성이 줄어들고 시스템의 나머지 부분과 독립적으로 코드를 검증 할 수 있으며 경우에 따라 자체 검증 테스트를 수행 할 필요가 있습니다. Test Double은 이러한 객체에 사용되는 일반적인 용어입니다.

이러한 독립된 이 오브젝트를 테스트 용으로 분리하는  일반적으로 테스트 더블이라 불린다.

https://www.slideshare.net/youngeunchoi12/effective-unit-testing-ch3

https://blog.pragmatists.com/test-doubles-fakes-mocks-and-stubs-1a7491dfa3da

테스트 더블은 테스트 할 대상인 오브젝트(클래스)를  영화에서 나오는 스턴트 대역(Stunt double)에서 유례된 말이다. 테스트시에 실제 객체를 대신 할 수 있는 객체를 의미 합니다. (여기서 ‘더블’이란 말의 유래는 영화에서 스턴트 대역[stunt double]을 생각하시면 될 듯 하네요.)

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

Doummy object는 전달되지만 실제로 사용되지는 않습니다. 일반적으로 매개 변수 목록을 채우기 위해 사용됩니다.

Fake object는 실제로 실제로 구현되어 있지만 일반적으로 제작에 적합하지 않은 단축키를 사용합니다 ( InMemoryTestDatabase 가 좋은 예입니다).

Test Stub 은 테스트 중에 작성된 호출에 대한 미리 준비된 답변을 제공합니다. 일반적으로 테스트를 위해 프로그래밍 된 내용 외에는 응답하지 않습니다.

Test Spy 은 그들이 불리는 방법에 따라 정보를 기록하는 스텁입니다. 한 가지 형태는 전송 된 메시지의 수를 기록하는 전자 메일 서비스 일 수 있습니다.

Mocks 는 그들이 기대하는 호출의 명세를 형성하는 기대치로 미리 프로그램되어있다. 기대하지 않은 전화를 받았을 때 예외를 throw 할 수 있으며 확인 중에 전화를 통해 예상했던 모든 전화를 받았는지 확인할 수 있습니다.


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



https://martinfowler.com/bliki/TestDouble.html

https://github.com/testdouble/contributing-tests/wiki/Test-Double

https://lostechies.com/derekgreer/2011/05/15/effective-tests-test-doubles/

https://laurentkempe.com/2010/07/17/Unit-Test-using-test-doubles-aka-Mock-Stub-Fake-Dummy/

[https://en.wikipedia.org/wiki/Test_double](https://en.wikipedia.org/wiki/Test_double)

[https://stackoverflow.com/questions/12827580/mocking-vs-spying-in-mocking-frameworks](https://stackoverflow.com/questions/12827580/mocking-vs-spying-in-mocking-frameworks)

[https://eminentstar.github.io/2017/07/24/about-mock-test.html](https://eminentstar.github.io/2017/07/24/about-mock-test.html)

[http://www.jpstory.net/2013/07/26/know-your-test-doubles/](http://www.jpstory.net/2013/07/26/know-your-test-doubles/)