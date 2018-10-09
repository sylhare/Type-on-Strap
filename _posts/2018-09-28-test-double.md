---
layout: post
title: "TDD : Test Double"
tags: [TDD, TestDouble, Mock, Stub, Fake]
#display: "false"
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

 테스트 더블이라는 단어의 기원은 Stunt Double에서 비롯됐다.
 
스턴트 더블은 영화에서 사용되는 용어인데 이 용어를 테스트에서도 착안하여 사용하고 있는 이유는 스턴트 더블의 배경을 알면 이해하기 쉽다.

_a trained professional who stands in for an actor in order to perform dangerous or physically demanding stunts._

 영화에서 건물 혹은 차량이 폭파하거나 차량 추격전 같은 위험한 장면이 있다면 배우는 영화를 위해 이러한 상황을 극복하고 연기를 해야 하지만 연기를 하다가 상처를 입거나 심한 경우 영화에 차질이 생겨 영화 자체가 망하는 위험이 생길 수 있다.

<img src="/md/img/test-double/stunt-double.jpg">
<em>Actor with stunt double</em>

물론 위험한 장면이 필요하다면 감독은 해당 장면을 찍어야 하지만 따르는 위험을 측정하기 어렵고 두렵다. 또한, 영화 투자자들도 위험을 좋아하지 않으므로, 이에 감독은 배우에게 위험한 장면을 연기하라고 요청하는 대신 배우처럼 보이는 다른 사람(스턴트 더블)에게 해당 연기를 시킨다. 이처럼 영화에선 여러 위험요소를 최소화하기 위해 스턴트 더블을 사용한다.

영화에서 스턴트 더블이 있다면 테스트엔 테스트 더블이 있다. 테스트 더블은 실제 객체와 같은 동작을 하는 객체를 뜻한다.

---

### Independence

프로그램도 마찬가지다. 물론 때에 따라 실제 객체를 대상으로 테스트를 수행하는 때도 있지만, 일반적으로 테스트 시에 해당한 객체를 직접 테스트하면 위험요소가 따른다. 
그 이유인즉슨 테스트 대상이 될 객체를 직접 사용하게 되면 실제 객체와 관계가 형성된 객체 간의 상호관계에서 발생하는 부작용을 일으키기 때문이다.

<img src="/md/img/test-double/test-real-object.png">
<em>Class Dependency Relationship</em>

다음 그림과 같이 테스트 시 실제 객체를 사용하면 이와 관계를 맺은 다른 객체들에 의해 값이 변경되거나 해당 객체만 테스트하기 어려운 부작용이 생긴다. 이러한 부작용들 때문에 해당 객체를 독립적인 객체로 만드는 코드를 작성하는 경우가 생긴다.
이러한 행위는 결과적으로 추가적인 코드로 인해 테스트 코드가 복잡해진다. 

<img src="/md/img/test-double/test-double-object.png">
<em>Class Independency Relationship of Test Double</em>

테스트 더블은 다음 그림과 같이 관계를 맺은 객체를 배제해 독립성을 갖게 해주고 실체 객체와 같은 행동을 하는 객체를 생성한다.
이처럼 독립된 객체를 테스트로 사용하는 이유는 테스트 시에 발생하는 예측 불가능한 위험요소를 최소화할 방법이기 때문이다. 또한, 테스트 코드의 복잡성이 줄이고 시스템의 나머지 부분과 독립적으로 코드를 검증할 수 있게 도와준다. 이외에 특수한 상황을 테스트한다거나 감춰진 정보를 얻기 위해 사용하기도 한다.

테스트 더블에는 `Test stub`, `Mock object`, `Test spy`, `Fake object`, `Dummy object`가 있다. 각각의 용도가 다르므로 테스트 객체의 목적에 따라 구분하고 사용해야 한다.

> [Types of test doubles](https://en.wikipedia.org/wiki/Test_double) <br/>
> - Test stub <br/>
> - Mock object <br/>
> - Test spy <br/>
> - Fake object <br/>
> - Dummy object <br/>

---

### Stub

_used for providing the tested code with "indirect input"_

Test Stub은 사전에 정의된 데이터를 보유하고 특정 객체 호출할 때 값을 정의해둔 데이터로 대체하는 테스트 더블이다. 즉 Stub은 로직이 없고 단지 원하는 값을 반환한다. 이 때문에 테스트 시에 '이 객체는 무조건 이 값을 반환한다.'라고 가정할 경우 사용한다. 또한, 실제 객체의 데이터를 사용할 때 위험이 있거나 객체를 포함할 수 없는 경우, 또는 원하지 않을 때 사용한다.

<img src="/md/img/test-double/stub.png">
<em>Test Stub</em>

가장 단순한 예로는 메소드 호출하여 데이터베이스에서 데이터를 가져와야 하는 경우다. 이 경우 실제 데이터베이스와 통신하는 객체 대신 Stub을 통해 반환할 데이터를 정의한다.

``` java
public class BankService {
    private final BankFactor bankFactor;
    
    public BankService(BankFactor bankFactor) {
        this.bankFactor = bankFactor;
    }
    
    Double getAvgWage(UserDAO userDao) {
        return calAvgWage(bankFactor.selectUserAmt(userDao));
    }
}
```

다음 코드는 데이터베이스에 접근하여 반환된 값을 통해 사용자의 평균임금을 계산하는 로직이다.

해당 서비스 클래스를 테스트하기 위해선 이와 . 단지 테스트를 진행하려 했을 뿐인데 배꼽이 배보다 큰 상황이다.
개발자는 서비스 클래스만 테스트하면 될 일을 팩토리 클래스까지 검증해야 할 상황에 직면한다. 이 모든 시작은 데이터베이스에서 데이터를 얻기 위해 시작됐다. 즉 테스트 시 `selectUserAvgWage` 메소드를 직접 테스트하기엔 위험이 따른다.

Stub는  데이터 베이스를 접근하여 사용하기 전에 대신 해당 데이터를 Stub으로 사전 구성한다. 데이터는 해당 알고리즘을 고려하고 테스트하기에 충분한 데이터로 정의한다.

``` java
public class BankServiceTest {
    private UserDAO    userDao;
    private BankFactor bankFactor;
	private StubData   stubData;
	
    @Before
    public void setUp() throws Exception {
        bankFactor = mock(BankFactor.class);
        userDao  = new UserDAO();
        StubData = new StubData();
    }
    
    @Test
    public void getAvgWageTest() {
        when(bankFactor.selectUserAmt(userDao))
        	.thenReturn(stubData.avgWage()); //stubbing bankFactor
        	
        double avgWage = new BankService(bankFactor).getAvgWage(userDao);
        assertThat(avgWage).isEqualTo(1000000.0);
    }
}
```



결과를 반환하고 시스템의 상태를 변경하지 않는 메소드를 조회 라고 합니다. 방법 averageGrades 학생 성적의 평균 반환 좋은 예입니다.
Double averageGrades(Student student);

값을 반환하고 부작용이 없습니다. 학생들의 채점 예제에서 보았 듯이이 유형의 메소드를 테스트하기 위해 스텁을 사용합니다. 우리는 실제 기능을 대체하여 메소드가 작업을 수행하는 데 필요한 값을 제공합니다. 그런 다음, 메소드에 의해 리턴 된 값은 표명에 사용될 수 있습니다.

명령 이라는 또 다른 카테고리가 있습니다 . 이것은 메소드가 시스템 상태를 변경하는 일부 조치를 수행하지만 리턴 값을 기대하지 않습니다.
void sendReminderEmail(Student student);

좋은 방법은 객체의 메서드를 두 개의 분리 된 범주로 나누는 것입니다. 
이 사례의 이름은 다음과 같습니다 : Bertrand Meyer의 "객체 지향 소프트웨어 구성" 에서 명령 쿼리 분리 .

Query 유형 메소드를 테스트하기 위해 우리는 메소드의 리턴 값을 검증 할 수 있기 때문에 스텁 사용을 선호해야합니다. 하지만 전자 메일을 보내는 방법과 같은 메서드의 명령 유형은 어떻습니까? 값을 반환하지 않을 때 테스트하는 방법? 대답은 모의 (Mock)입니다. 우리가 다룰 마지막 유형의 테스트 더미입니다.

=======
>>>>>>> branch 'master' of https://github.com/gmun/gmun.github.io

<img src="/md/img/test-double/stub.png">
<em>Test Stub</em>

가장 단순한 예로는 메소드 호출하여 데이터베이스에서 데이터를 가져와야 하는 경우다. 이 경우 실제 데이터베이스와 통신하는 객체 대신 Stub을 통해 반환할 데이터를 정의한다.

``` java
public class BankService {
    private final BankFactor bankFactor;
    
    public BankService(BankFactor bankFactor) {
        this.bankFactor = bankFactor;
    }
    
    Double getAvgWage(UserDAO userDao) {
        return calAvgWage(bankFactor.selectUserAmt(userDao));
    }
}
```

다음 코드는 데이터베이스에 접근하여 반환된 값을 통해 사용자의 평균임금을 계산하는 로직이다.

해당 서비스 클래스를 테스트하기 위해선 이와 . 단지 테스트를 진행하려 했을 뿐인데 배꼽이 배보다 큰 상황이다.
개발자는 서비스 클래스만 테스트하면 될 일을 팩토리 클래스까지 검증해야 할 상황에 직면한다. 이 모든 시작은 데이터베이스에서 데이터를 얻기 위해 시작됐다. 즉 테스트 시 `selectUserAvgWage` 메소드를 직접 테스트하기엔 위험이 따른다.

Stub는  데이터 베이스를 접근하여 사용하기 전에 대신 해당 데이터를 Stub으로 사전 구성한다. 데이터는 해당 알고리즘을 고려하고 테스트하기에 충분한 데이터로 정의한다.

``` java
public class BankServiceTest {
    private UserDAO    userDao;
    private BankFactor bankFactor;
	private StubData   stubData;
	
    @Before
    public void setUp() throws Exception {
        bankFactor = mock(BankFactor.class);
        userDao  = new UserDAO();
        StubData = new StubData();
    }
    
    @Test
    public void getAvgWageTest() {
        when(bankFactor.selectUserAmt(userDao))
        	.thenReturn(stubData.avgWage()); //stubbing bankFactor
        	
        double avgWage = new BankService(bankFactor).getAvgWage(userDao);
        assertThat(avgWage).isEqualTo(1000000.0);
    }
}
```



결과를 반환하고 시스템의 상태를 변경하지 않는 메소드를 조회 라고 합니다. 방법 averageGrades 학생 성적의 평균 반환 좋은 예입니다.
Double averageGrades(Student student);

값을 반환하고 부작용이 없습니다. 학생들의 채점 예제에서 보았 듯이이 유형의 메소드를 테스트하기 위해 스텁을 사용합니다. 우리는 실제 기능을 대체하여 메소드가 작업을 수행하는 데 필요한 값을 제공합니다. 그런 다음, 메소드에 의해 리턴 된 값은 표명에 사용될 수 있습니다.

명령 이라는 또 다른 카테고리가 있습니다 . 이것은 메소드가 시스템 상태를 변경하는 일부 조치를 수행하지만 리턴 값을 기대하지 않습니다.
void sendReminderEmail(Student student);

좋은 방법은 객체의 메서드를 두 개의 분리 된 범주로 나누는 것입니다. 
이 사례의 이름은 다음과 같습니다 : Bertrand Meyer의 "객체 지향 소프트웨어 구성" 에서 명령 쿼리 분리 .

Query 유형 메소드를 테스트하기 위해 우리는 메소드의 리턴 값을 검증 할 수 있기 때문에 스텁 사용을 선호해야합니다. 하지만 전자 메일을 보내는 방법과 같은 메서드의 명령 유형은 어떻습니까? 값을 반환하지 않을 때 테스트하는 방법? 대답은 모의 (Mock)입니다. 우리가 다룰 마지막 유형의 테스트 더미입니다.



---

### Mock

_used for verifying "indirect output" of the tested code, by first defining the expectations before the tested code is executed_

Mocks 는 그들이 기대하는 호출의 명세를 형성하는 기대치로 미리 프로그램되어있다. 기대하지 않은 전화를 받았을 때 예외를 throw 할 수 있으며 확인 중에 전화를 통해 예상했던 모든 전화를 받았는지 확인할 수 있습니다.

Mock은 수신 호출을 등록하는 객체입니다. 
테스트 주장에서 우리는 모의 (Mock)에서 모든 예상 된 행동이 수행되었다는 것을 검증 할 수있다.

생산 코드를 호출하고 싶지 않거나 손쉬운 검증 방법이없는 경우 의도 된 코드가 실행되었음을 나타 내기 위해 모의 (mock)를 사용합니다. 반환 값은 없으며 시스템 상태 변경을 확인하는 쉬운 방법이 없습니다. 전자 메일 전송 서비스를 호출하는 기능을 예로들 수 있습니다. 
우리는 테스트를 실행할 때마다 전자 메일을 보내지 않습니다. 더욱이, 테스트에서 올바른 이메일이 전송되었는지 확인하는 것은 쉽지 않습니다. 우리가 할 수있는 일은 테스트에서 수행 된 기능의 출력을 검증하는 것입니다. 다른 세계에서는 전자 메일 보내기 서비스가 호출되었는지 확인하십시오.

다음 예제에서 비슷한 경우가 나타납니다.

<img src="/md/img/test-double/mock.png">
<em>Mock object</em>



---

### Spy

_used for verifying "indirect output" of the tested code, by asserting the expectations afterwards, without having defined the expectations before the tested code is executed. It helps in recording information about the indirect object created_

Test Spy 은 그들이 불리는 방법에 따라 정보를 기록하는 스텁입니다. 한 가지 형태는 전송 된 메시지의 수를 기록하는 전자 메일 서비스 일 수 있습니다.

---

### Fake

_used as a simpler implementation, e.g. using an in-memory database in the tests instead of doing real database access_

Fake object는 실제로 실제로 구현되어 있지만 일반적으로 제작에 적합하지 않은 단축키를 사용합니다 ( InMemoryTestDatabase 가 좋은 예입니다).

가짜 (fake)는 작동하는 구현을 가지고 있지만 제작과는 다른 객체입니다. 보통 그들은 간단한 코드를 사용하고 생산 코드의 버전을 단순화했습니다.

이 바로 가기의 예는 데이터 액세스 개체 또는 저장소의 메모리 내 구현 일 수 있습니다. 이 가짜 구현은 데이터베이스에 관여하지 않지만 단순 콜렉션을 사용하여 데이터를 저장합니다. 이를 통해 데이터베이스를 시작하고 시간이 많이 소요되는 요청을 수행하지 않고도 서비스의 통합 테스트를 수행 할 수 있습니다.

<img src="/md/img/test-double/fake.png">
<em>Fake object</em>

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