---
layout: post
title: "TDD : Test Double"
tags: [TDD, TestDouble, Mock, Stub, Fake]
display: "false"
feature-img: "md/img/test-double/test-double-thumbnail.jpg"              
thumbnail: "md/img/test-double/test-double-thumbnail.jpg"
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

 테스트 더블의 기원은 Stunt Double에서 비롯됐다.
 
스턴트 더블은 영화에서 사용되는 용어인데 이 용어를 테스트에서도 착안하여 사용하고 있는 이유는 스턴트 더블의 배경을 알면 이해하기 쉽다.

_a trained professional who stands in for an actor in order to perform dangerous or physically demanding stunts._

 영화에서 건물 혹은 차량이 폭파하거나 차량 추격전 같은 위험한 장면이 있다면 배우는 영화를 위해 이러한 상황을 극복하고 연기를 해야 하지만 연기를 하다가 상처를 입거나 심한 경우 영화에 차질이 생겨 영화 자체가 망하는 위험이 생길 수 있다.

<img src="/md/img/test-double/stunt-double.jpg">
<em>Actor with stunt double</em>

물론 위험한 장면이 필요하다면 감독은 해당 장면을 찍어야 하지만 따르는 위험을 측정하기 어렵고 두렵다. 또한, 영화 투자자들도 위험을 좋아하지 않으므로, 이에 감독은 배우에게 위험한 장면을 연기하라고 요청하는 대신 배우처럼 보이는 다른 사람(스턴트 더블)에게 해당 연기를 시킨다. 이처럼 영화에선 여러 위험요소를 최소화하기 위해 스턴트 더블을 사용한다.

영화에서 스턴트 더블이 있다면 테스트엔 테스트 더블이 있다. 이처럼  테스트 더블은 실제 객체와 같은 동작을 하는 객체를 뜻한다.

---

### Independence

프로그램도 마찬가지다. 물론 때에 따라 실제 객체를 대상으로 테스트를 수행할 수 있지만, 일반적으로 테스트 시에 해당한 객체를 직접 테스트하면 위험요소가 따른다.

<img src="/md/img/test-double/test-real-object.png" style="max-height: 400px;">
<em>테스트 의존관계</em>

그 이유인즉슨 테스트 대상이 될 객체를 직접 사용하게 되면 실제 객체(Real Object)와 의존관계가 형성된 다른 객체 간의 상호관계에서 발생하는 부작용을 일으키기 때문에 테스트 시에 실체 객체를 대신 할 수 있는 독립된 객체(테스트 더블)를 생성하여 사용한다.

이처럼 테스트 더블을 사용하는 가장 큰 이유는 테스트 객체의 독립성을 보장받기 때문이다. 이는 테스트 시에 발생하는 위험요소를 최소화할 방법이다. 또한, 테스트 코드의 복잡성이 줄이고 시스템의 나머지 부분과 독립적으로 코드를 검증할 수 있게 도와준다.

테스트 더블에는 `Test stub`, `Mock object`, `Test spy`, `Fake object`, `Dummy object`가 있다.

> [Types of test doubles](https://en.wikipedia.org/wiki/Test_double) <br/>
> - Test stub <br/>
> - Mock object <br/>
> - Test spy <br/>
> - Fake object <br/>
> - Dummy object <br/>

테스트 더블의 종류는 각각의 용도가 다르므로 테스트 객체의 목적에 따라 구분하고 사용해야 한다.

---

### Stub

_used for providing the tested code with "indirect input"_

Test Stub 은 테스트 중에 작성된 호출에 대한 미리 준비된 답변을 제공한다. 일반적으로 테스트를 위해 프로그래밍 된 내용 외에는 응답하지 않는다.

Stub은 로직이 없고 단지 원하는 값을 반환합니다. 테스트시에 “이 객체는 무조건 이 값을 반환한다”고 가정할 경우 사용할 수 있습니다. Stub은 보통 작성하기 쉽지만 불필요한 boilerplate 코드를 줄이기 위해서 Mocking Framework을 이용하는게 편합니다.

---

### Mock

_used for verifying "indirect output" of the tested code, by first defining the expectations before the tested code is executed_

Mocks 는 그들이 기대하는 호출의 명세를 형성하는 기대치로 미리 프로그램되어있다. 기대하지 않은 전화를 받았을 때 예외를 throw 할 수 있으며 확인 중에 전화를 통해 예상했던 모든 전화를 받았는지 확인할 수 있습니다.


---

### Spy

_used for verifying "indirect output" of the tested code, by asserting the expectations afterwards, without having defined the expectations before the tested code is executed. It helps in recording information about the indirect object created_

Test Spy 은 그들이 불리는 방법에 따라 정보를 기록하는 스텁입니다. 한 가지 형태는 전송 된 메시지의 수를 기록하는 전자 메일 서비스 일 수 있습니다.

---

### Fake

_used as a simpler implementation, e.g. using an in-memory database in the tests instead of doing real database access_

Fake object는 실제로 실제로 구현되어 있지만 일반적으로 제작에 적합하지 않은 단축키를 사용합니다 ( InMemoryTestDatabase 가 좋은 예입니다).

---

### Dummy

_used when a parameter is needed for the tested method but without actually needing to use the parameter_

Doummy object는 전달되지만 실제로 사용되지는 않습니다. 일반적으로 매개 변수 목록을 채우기 위해 사용됩니다.

---

### 마무리

http://engineering.pivotal.io/post/the-test-double-rule-of-thumb/

https://www.slideshare.net/youngeunchoi12/effective-unit-testing-ch3

https://blog.pragmatists.com/test-doubles-fakes-mocks-and-stubs-1a7491dfa3da

http://xunitpatterns.com/Test%20Double.html


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