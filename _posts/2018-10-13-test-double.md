---
layout: post
title: "TDD : Test Double"
tags: [TDD, TestDouble, Mock, Stub, Fake]
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

프로그램도 마찬가지다. 물론 때에 따라 실제 객체를 대상으로 테스트를 수행하는 때도 있지만, 일반적으로 테스트 시에 실제 객체를 직접 테스트하면 실제 객체와 관계가 형성된 객체 간의 상호관계에서 발생하는 예측 불가능한 위험이 발생할 수 있다. 

<img src="/md/img/test-double/test-real-object.png">
<em>Class Dependency Relationship</em>

TestObject에 대한 테스트 코드를 작성한다고 가정하자. 여러 테스트 케이스를 적용하여 테스트 진행하던 중 에러가 발생하였다. 이때 Real Object는 불완전한 상태라고 가정한다면 개발자는 Real Object와 관계가 형성된 객체들부터 검증한 다음에 테스트를 진행해야 TestObject를 검증할 수 있다. 

또 다른 예로는 테스트 시 실제 객체와 관계를 맺은 다른 객체들에 의해 값이 변경되는 부작용이 발생할 수 있다. 이러한 예측 불가능한 위험 때문에 해당 객체를 독립적인 객체로 대체하여 테스트를 진행한다. 이러한 행위는 결과적으로 개발자는 많은 상황을 고려해야 하며 추가적인 코드로 인해 테스트 코드가 복잡해진다. 

<img src="/md/img/test-double/test-double-object.png">
<em>Class Independency Relationship of Test Double</em>

테스트 더블은 다음 그림과 같이 관계를 맺은 객체를 배제해 독립성을 갖게 해주고 실체 객체와 같은 행동을 하는 객체를 생성한다.

이처럼 독립된 객체를 테스트로 사용하는 이유는 테스트 시에 발생하는 예측 불가능한 위험요소를 최소화할 방법이기 때문이다. 또한, 테스트 코드의 복잡성이 줄이고 시스템의 나머지 부분과 독립적으로 코드를 검증할 수 있게 도와준다. 이외에 특수한 상황을 테스트한다거나 감춰진 정보를 얻기 위해 사용한다.

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

Test Stub은 사전에 정의된 데이터를 보유하고 특정 객체 호출할 때 값을 정의해둔 데이터로 대체하는 테스트 더블이다. 즉 Stub은 로직이 없고 단지 원하는 값을 반환한다. 이 때문에 테스트 시에 '이 객체는 무조건 이 데이터를 반환한다.'라고 가정할 경우 사용한다. 또한, 실제 객체의 데이터를 사용할 때 위험이 있거나 객체를 포함할 수 없는 경우, 또는 원하지 않을 때 사용한다.

<img src="/md/img/test-double/stub.png">
<em>Test Stub</em>

가장 단순한 예로는 메소드 호출하여 데이터베이스에서 데이터를 가져와야 하는 경우다. 이 경우 실제 데이터베이스와 통신하는 객체 대신 Stub을 통해 반환할 데이터를 정의한다.

``` java
public class BankService {
    private final BankFactor bankFactor;
    
    public BankService(BankFactor bankFactor){
        this.bankFactor = bankFactor;
    }
    
    public double getAvgWage(HashMap<String, Object> params){
        return calAvgWage(bankFactor.selectUserAmt(params));
    }
    ...
}

```

`selectUserAmt` 메소드를 호출하여 데이터베이스가 반환된 데이터를 사용하여 평균임금을 계산한다고 가정하자.

BankService의 테스트 코드를 작성 시 실제 데이터베이스를 접속해야 하고 BankFactor가 불완전한 상태라면 BankFactor 클래스부터 검증한 뒤에 테스트 코드를 작성해야 한다. 단지 서비스 클래스를 검증하려 했을 뿐인데 배꼽이 배보다 큰 상황이다. 이 모든 시작은 서비스 클래스와 팩토리 클래스가 의존 관계이기 때문이다.

독립적인 서비스 클래스를 테스트하기 위해 팩토리 클래스가 데이터베이스에 접근하여 반환한 데이터를 Stub으로 사전에 정의하여 해당 클래스를 검증해보자.

``` java
public class BankServiceTest {
    private BankFactor bankFactor;
	
    @Before
    public void setUp() throws Exception {
        bankFactor = mock(BankFactor.class);
    }
    
    @Test
    public void getAvgWage(){
        when(bankFactor.selectUserAmt(new HashMap<String, Object>()))
        	.thenReturn(new StubDatas().userAmts()); //Stubbing bankFactor
        	
        double avgWage = new BankService(bankFactor).getAvgWage(new HashMap<String, Object>());
        assertThat(avgWage).isEqualTo(1000000.0);
    }
}
```

다음 테스트 코드처럼 `bankFactor.selectUserAmt`의 반환 값 대신 평균임금 계산 알고리즘을 고려하여 사전에 정의한 메소드의 결괏값을 사용하여 검증한다.

`.thenReturn(new StubDatas().userAmts());`

실제 객체를 대신하여 StubDatas()를 사용함으로써 테스트를 위해 BankFactor 클래스를 수정하거나 코드를 추가하는 별도의 작업을 할 필요가 없고 반환되는 값에 대한 부작용을 고려하지 않아도 된다.
이처럼 로직을 테스트하는 데 필요한 값을 Stub으로부터 받아 검증한다.

`void sendResultAvgWageEmail(UserDAO user, String content);`

다음과 같이 계산된 평균임금 값을 사용자의 메일로 발송하는 메소드가 있다고 가정하자. 반환된 값이 있다면 Stub을 통해 검증하면 되지만 반환하는 값이 없다면 어떻게 검증해야 할까? 이런 경우 Mock을 통해 검증하면 된다.

---

### Mock

_used for verifying "indirect output" of the tested code, by first defining the expectations before the tested code is executed_

Mock은 호출에 대한 기대하는 실행 결과를 사전에 정의한 객체다. 이 객체를 통해 기대하지 않은 결과를 예외 처리할 수 있으며 예상했던 모든 결과를 확인할 수도 있다. 
일반적으로 실제 코드를 호출하고 싶지 않거나 손쉬운 검증 방법이 없는 경우 의도된 코드가 실행되었음을 나타내기 위해 사용한다. 즉 Mock은 동작에 대한 검증으로 반환 값은 없다. 동작에 대한 검증은 테스트할 수 있지만, 동작하는 그 자체를 검증하는 것은 어렵다.

<img src="/md/img/test-double/mock.png">
<em>Mock Object</em>

단편적인 예로는 메일 서비스 기능을 들 수 있다. BankService에서 사용자의 평균임금을 계산한 값을 사용자의 메일로 전송해주는 기능을 추가되었다고 가정하자. 테스트 코드를 실행할 때마다 결괏값이 메일로 전송하지 않는다. 더구나 메일이 올바르게 전송되었는지 확인하기는 쉽지 않다. 개발자가 할 수 있는 일은 테스트에서 수행된 기능의 출력을 검증하는 것이다.

``` java
public class BankService {
    private final BankFactor bankFactor;
    private final MailService mailService;
    
    public BankService(BankFactor bankFactor, MailServiceImple mailServiceImple){
    	this.bankFactor = bankFactor;
    	this.mailService = mailServiceImple;
    }
    
    public double getAvgWage(HashMap<String, Object> params){
    	List<UserDAO> users = bankFactor.selectUserAmt(params);
    	double avgWage = calAvgWage(users);
    	
    	mailService.sendResultAvgWageEmail(users, Double.toString(avgWage));
        return avgWage;
    }
}
```

다음 코드를 보면 BankService에 사용자의 요구대로 결과 값을 평균임금 값을 메일로 전송해주는 기능을 추가했다.

`mailService.sendResultAvgWageEmail(users, Double.toString(avgWage));`

테스트 시 테스트 결과 값을 메일로 전송하지 말아야 한다. 메일 기능을 담고 있는 MailService 클래스를 Mock 객체로 배치하여 이를 해결한다.

``` java
public class BankServiceTest {
    private BankFactor bankFactor;
    private MailServiceImple mailServiceMock;
    
    @Before
    public void setUp() throws Exception {
        bankFactor = mock(BankFactor.class);
        mailServiceMock = mock(MailServiceImple.class);
    }
    
    @Test
    public void getAvgWage() throws Exception{
    	List<UserDAO> data = new StubDatas().userDaoList();
    	
        when(bankFactor.selectUserAmt(new HashMap<String, Object>()))
        	.thenReturn(data); //Stubbing bankFactor
        	
        double avgWage = new BankService(bankFactor, mailServiceMock).getAvgWage(new HashMap<String, Object>());
        
        verify(mailServiceMock).sendResultAvgWageEmail(data, Double.toString(avgWage));
        assertThat(avgWage).isEqualTo(1000000.0);
    }
}
```

`getAvgWage` 메소드 실행 후 `sendResultAvgWageEmail` 메일 메소드가 실행되었는지 verify 메소드를 통해 확인했다. 하지만 사용자에게 제대로 메일이 전달되었는지 해당 기능이 예외를 반환하는지에 대한 검증은 확인할 수 없다. 이것이 MailService 관점에서 테스트할 수 있는 전부다.

`verify(mailServiceMock).sendResultAvgWageEmail(data, Double.toString(avgWage));`

여기서 '실제 사용자에게 제대로 메일이 전달되었는지 알 수는 없을까?'라는 의문이 든다. 그에 대한 답은 알 수 없다. 그러나 현시점에서 그 의문에 대해 신경 쓰지 않아도 된다. 이는 BankService의 책임이 아닌 MailService 책임이기 때문이다.

---

### Spy

_used for verifying "indirect output" of the tested code, by asserting the expectations afterwards, without having defined the expectations before the tested code is executed. It helps in recording information about the indirect object created_

Test Spy는 실제 객체의 메소드를 호출하고 반환 값이 있으면 해당 반환값도 반환해준다. 

<img src="/md/img/test-double/spy.png">
<em>Test Spy</em>

일반적으로 Spy는 시스템이 메소드를 호출했는지 확인하고 싶을 때 사용한다. 예를들어 호출 횟수를 계산하거나 전달 된 인수를 기록하는 것과 같은 모든 종류의 것을 기록하는 목적인 경우 유용하다.

``` java
@Test
public void userBankCountTest() throws Exception {
	List<UserDAO> bankList = bankFactor.selectFindByBankName("KR은행"); 
	List<UserDAO> spy = spy(bankList);
	
	when(spy.size()).thenReturn(5); //stubbing list size
	
	spy.add(new UserDAO(1, "A"));
	spy.add(new UserDAO(2, "B"));

	System.out.println(spy.get(0)); //A
	System.out.println(spy.size()); //5

	verify(spy).add(new UserDAO(1, "A"));
	verify(spy).add(new UserDAO(2, "B"));

	
	when(spy.get(100))
		.thenReturn(new UserDAO(1, "A")); // IndexOutOfBoundsException 
}
```

해당 메소드는 사용자가 등록한 은행의 개수를 반환한다. Mockito 프레임워크에서 제공하는 spy 메소드를 사용하여 반환한 객체를 Spy했다. 당연히 spy된 객체를 Stub을 할 수 있다.

`doReturn(new UserDAO(1, "A").when(spy).get(100);`

단 Spy를 활용하여 Stub을 하면 실제 인스턴스의 메소드를 호출하기 때문에 종종 예기치 못한 예외가 발생한다. 위의 when(spy.get(100))에서는 진짜 인스턴스의 메소드를 호출하기 때문에 IndexOutOfBountException이 발생하게 된다. 이 경우 Mockito.doReturn를 사용해서 문제를 회피할 수 있다.

``` java
public class MailServiceImpleTest {
	private SpyFileIO spy;
	
	@Before
	public void setUp(){
		spy = spy(SpyFileIO.class);
	}
	
	@Test
	public void sendTest(){
		MailServiceImple mailSvc = new MailServiceImple(spy);
		
		mailSvc.send("gmun0929@gmail.com", "제목", "내용", null);
		verify(mailSvc).send(null, null, null, null);
		assertEquals(1, spy.callCount);
	}
	
	private class SpyFileIO implements FileIO{
		public int callCount = 0;
		
		@Override
		public StringBuilder read(String filePath){
			this.callCount++;
	        return null;
		}
	}
}
```

다음과 같이 Spy를 활용하여 테스트 시 특정 메소드가 호출된 총횟수를 검사할 수도 있다.


---

### Fake

_used as a simpler implementation, e.g. using an in-memory database in the tests instead of doing real database access_

Fake는 실제 데이터베이스의 데이터에 접근하는 객체, 즉 실제 DAO를 대신하여 테스트를 수행할 객체를 의미한다.

<img src="/md/img/test-double/fake.png">
<em>Fake Object</em>

Fake는 실제 데이터베이스가 응답한 데이터의 축소판이라고 생각하면 된다. Fake 구현은 실제 데이터베이스에 관여하지 않고 단순 Collection을 사용하여 테스트 시 필요한 데이터를 저장한다. 이를 통해 데이터베이스의 응답, 요청 시간이 많이 소요되는 서비스 통합 테스트를 비교적 빠르게 검증할 수 있다.

``` java
public class FakeBankRepository {
	private UserDAO user = new UserDAO();
	private BankDAO bank = new BankDAO();
	private Map<UserDAO, BankDAO> userAmts = new HashMap<>();
    
	public FakeBankRepository(){
		this.user.setId(100L);
		this.bank.setBankName("KR은행");
		this.userAmts.put(user, bank);
	}
     
	String getUserBankName(UserDAO user){
		return userAmts.get(user).getBankName();
	}
}
```

실제 DAO를 대체할 FakeBankRepository 객체를 생성한 뒤 테스트에 필요한 정보를 Collection에 담아 테스트 시 필요한 데이터 Collection을 통해 받는다.

이때 Fake와 실제 DAO의 테스트에 관한 결과는 같을 수 있지만, 엄밀히 따지면 다르다. 이 문제점은 개발자들이 흔히 실수하는 부분인데 객관적이면서 축소된 데이터를 가지고 기능을 검증한 뒤 기능이 완벽하다고 단정을 짓는 것이다.

이러한 문제점을 대체할 Fake 구현 방식은 여러 방법이 있겠지만, 일반적으로 인메모리 데이터베이스의 DAO를 Fake로 사용한다.

``` java
@H2DB
public interface BankRepository extends JpaRepository<BankDAO, Long>{
	public List<BankDAO> findByIdIn(List<Long> ids);
}
```

스프링의 AOP를 사용하여 인메모리 데이터베이스로 접속하게 설정했다. 이 때문에 클래스 코드를 변경하지 않고 DAO는 인메모리 데이터베이스의 데이터에 영향을 받을 수 있게 됐다. 이때 해당 인메모리의 데이터는 테스트에 필요한 데이터만 삽입하여 사용한다.

하지만 인메모리 데이터베이스 또한 테스트의 정확성을 보장하진 않는다. 인메모리를 사용하는데 속도를 개선할 순 있지만, 인메모리를 사용한 테스트 결과에 대한 정확도는 실제 데이터를 검증한 정확도에 비해 떨어질 수밖에 없다. 때문에 테스트 시 Fake의 구현 방식을 권장하지 않고 필요하다면 실제 데이터베이스의 데이터를 통해 기능을 검증한다.

물론 Fake가 그림의 떡은 아니다. 에자일 개발 방법론에서 [스파이크](https://en.wikipedia.org/wiki/Spike_(software_development))나 사용자에게 빠른 피드백을 받기 위한  [프로토 타이핑](https://en.wikipedia.org/wiki/Prototype)에 활용한다면 유용하다. 즉 실제 데이터베이스 설계에 관한 결정을 미루고 인메모리 데이터베이스를 통해 시스템을 구현하고 실행할 수 있다는 장점을 잘 활용해야 한다.

---

### Dummy

_used when a parameter is needed for the tested method but without actually needing to use the parameter_

Dummy는 미국의 대표적인 코믹 영화 Dumb and Dumber(덤 앤 더머)와 연관하면 이해하기 쉽다.

영화 제목의 Dumb과 Dumber은 바보의 대명사처럼 모자란 사람을 부를 때 쓰였다. 이와 마찬가지로 Dummy라는 이름에서 알 수 있듯이 Dummy는 매우 바보 같은 객체다. 일반적으로 Dummy object는 해당 객체가 어떻게 사용되는지 상관없이 컴파일과 런타임 실행을 만족 시키기기 위해  객체를 전달할 때 사용한다. 즉 Dummy object는 일반적으로 매개 변수 목록을 채우기 위해 사용한다.

<img src="/md/img/test-double/dummy.png">
<em>Dummy Object</em>

 예를 들어 테스트 코드에 어느 한 객체가 매개변수가 있는 생성자를 포함하고 있다고 가정하자. 이때 매개 변수를 주입해야 하지만 해당 매개 변수는 테스트 시 해당 매개 변수를 사용하지 않는다면?
 
 이 경우에 Dummy를 사용하여 해결해보자.

``` java
public class MailServiceImple implements MailService{
	private final FileIO fileIO;
    
	public MailServiceImple(FileIO fileIO){
		this.fileIO = fileIO;
	}
    
	@Override
	public void send(String toEmail, String subject, String content, File attachFile){
		...
	}
    
	...
}
```

먼저 테스트 코드를 작성할 MailServiceImple를 보면 생성자를 통해  FileIO와 의존관계를 맺어주고 있다. 주입 받은 FileIO는 테스트 시 사용하지 않으리라고 간주하고 해당 객체를 Dummy하여 테스트 코드를 작성해보자.

```java
public class MailServiceImpleTest {
	private DummyFileIO dummyFileIO;

	@Before
	public void setUp(){
		dummyFileIO = mock(DummyFileIO.class);
	}
	
	@Test
	public void sendTest(){
		MailServiceImple mailSvc = new MailServiceImple(dummyFileIO);
        
		mailSvc.send("gmun0929@gmail.com", "제목", "내용", null);
		verify(mailSvc).send(null, null, null, null);
	}
    
	private class DummyFileIO implements FileIO{
		@Override
		public StringBuilder read(String filePath){
			throw new RuntimeException("Not expected to be called");
		}
	}
}
```
 
 눈치를 챘는지 모르겠지만 다른 테스트 더블의 구현 방식과는 다르게 Dummy인 경우 내부 클래스로 작성했다. 일반적으로 Dummy는 테스트를 통해 변경되지 않기 때문에 내부 클래스를 만들고 모든 테스트에 재사용하는 것이 더 적합하기 때문이다.

`MailServiceImple mailSvc = new MailServiceImple(null)`

 또한, 테스트 시 해당 인수 값이 무관하다면 null로 대체해도 무관하다.

`dummyFileIO = mock(DummyFileIO.class);`

파일 입출력 기능에서 발생하는 예외를 방지하기 위해 해당 인수 값을 mock 객체를 사용하여 독립적인 메일 발송 테스트를 진행했다.

---

### 마무리

<img src="/md/img/test-double/test-double-relationship.png">
<em>Goal of Test double use</em>

테스트 더블은 실제 객체와 관계를 맺은 객체들을 테스트용 객체로 대체하여 독립적인 테스트를 가능할 수 있게 만들어 주는 목적이 있다. 앞서 설명했듯이 테스트 더블에는 총 다섯 가지의 종류가 있고 이를 요약함으로써 글을 마친다.

`Test Stub`은 로직이 없고 사전에 정의한 데이터를 반환한다.

`Mock Object`는 객체를 동작이 없고 반환 값이 없는 상태로 만든다.

`Test Spy`는 실제 객체와 같은 동작을 한다.

`Fake Object`는 실제 DAO를 대체할 객체이다.

`Dummy Object`는 동작하지 않고 매개 변수 목록을 채워주기 위해 사용한다. 또한, 컴파일을 가능하게 하기 위한 것일 뿐이며 테스트에 포함되지 않는다.

---

### 참고

[xUnit - TestDouble](http://xunitpatterns.com/Test%20Double.html)<br/>
[Martin Fowler - TestDouble](https://martinfowler.com/bliki/TestDouble.html)<br/>
[Martin Fowler - Mocks Aren't Stubs](https://martinfowler.com/articles/mocksArentStubs.html)<br/>
[Martin Fowler - CommandQuerySeparation](https://martinfowler.com/bliki/CommandQuerySeparation.html)<br/>
