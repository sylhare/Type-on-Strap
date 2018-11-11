---
layout: post
title: "단위 테스트 코드를 작성 해보자"
tags: [TDD, TDDCycle, RGRCycle, JUnit, JUnit4, UnitTest]
categories: [Test]
subtitle: "RGR 단계를 인식하며 단위 테스트 코드를 개발하기"
feature-img: "md/img/thumbnail/tdd-cycle-practice.png"    
thumbnail: "md/img/thumbnail/tdd-cycle-practice.png"
excerpt_separator: <!--more-->
sitemap:
changefreq: daily
priority: 1.0
---

<!--more-->

# 단위 테스트 코드 입문기

---

본 포스팅은 단위 테스트 코드 개발에 처음 입문한 사람들을 대상으로한다.

이 글에선 개념적인 부분보다는 JUnit을 통해 단위 테스트를 작성하는 방법과 테스트 검증 시 일어나는 TDD-Cycle의 RGR 단계를 인지하여 개발하는 데에 초점을 둘 예정이다.

TDD-Cycle이 왜, 언제 발생되는지 모른다거나 자세한 개념적인 부분이 필요하다면 이전에 작성한 [TestFramework : JUnit](https://gmun.github.io/test/2018/11/04/junit.html)과 [TDD : Test-Driven Development](https://gmun.github.io/test/2018/08/24/test-driven-development.html)를 같이 참고한다면 본 포스팅을 읽는 데 도움이 될 것 같다.

### Download

JUnit 버전에 따라 요구하는 Java 버전이 다르지만 본 포스팅에선 JUnit4를 사용하기 때문에 Java7이면  충분하다. 

현재 JUnit4의 최종 버전은 4.12로 `junit.jar`, `hamcrest-core.jar` 두 가지 라이브러리를 필요로하다. [JUnit4.12 다운로드 홈페이지](https://github.com/junit-team/junit4/wiki/Download-and-Install)에선 JAR, Maven, Gradle의 설치방법을 가이드 해주고 있다. 현재 자신의 프로젝트가 라이브러리를 관리하는 방식에 따라 JUnit4 라이브러리를 추가하자.

>JUnit5는 최소 Java8 이상의 자바 버전을 요구한다. 이전 버전에서 지원하지 않았던 람다식, 스트림을 테스트할 수 있도록 지원하고 테스트 필터링, 확장 등 더 다양한 테스트 환경을 제공하고 있다.

### Test SourceFolder

프로젝트에 JUnit을 추가해줬다면 테스트 클래스를 관리할 소스 폴더를 만들어보자.

<img src="/md/img/unit-test-practice/test-source-folder.png" height="400px">
<em>Test SourceFolder</em>

소스 폴더 이름은 test로 만들었고 앞으로 작성된 테스트 클래스는 이 폴더에 관리할 예정이다. 이는 실제 코드와 테스트 코드를 따로 분리하여 관리함으로써 자연스레 기능적 응집도(Functional Cohesion)와 논리적 응집도(Logical Cohesion)를 높아지는 효과가 있다. 이 행위는 테스트 클래스를 관리를 쉽게 하고자 하는 목적이다.

### 주제 선정

 단위 테스트를 연습할 주제를 선정해보자. 가장 처음 [TDD](https://gmun.github.io/test/2018/08/24/test-driven-development.html)를 연습할 때는 유틸성 기능 또는 알고리즘 문제를 통해 연습하면 좋다.
 
> 문제 1) 더하기, 빼기, 나누기, 곱셈할 수 있는 문자열 계산기 만들기

다음의 문제는 유틸성 문제 중 가장 기초적인 `문자열 계산기 만들기`이다. 이를 통해 단위 테스트를 연습해보려 한다.

### TestCase를 작성하자

 기능의 난도가 높든 낮든 무턱대고 코드부터 짜는 습관은 매우 안 좋은 습관이다. 물론 테스트 케이스를 작성하는 것은 귀찮고 수고스럽지만, 이 단계는 매우 중요하다.
 
  일반적으로 주제의 목표에 도달하기 위해, 요구사항을 충분히 분석하고 테스트 케이스를 추출한다. 이러한 일련의 과정에 개발하기 전 기능의 전반적인 흐름을 알 수 있고 기존의 설계를 개선하여 더 나은 설계를 도출할 수도 있다.

<img src="/md/img/unit-test-practice/todo-list.png">
<em>Todo-List</em>
 
  테스트 케이스를 작성하는 이유는 이정표 역할을 하기 때문이다. 개발자는 테스트 케이스를 통해 개발의 진척도를 직관적으로 알 수 있고 목록을 하나씩 지울 때마다 소소한 성취감을 느낄 수 있다. 또한, 정리된 todo-list 덕분에 앞만 보고 개발할 수 있게 된다.

좋은 테스트 케이스를 작성하려면 기능의 목표를 설정하는 것이 중요하다. 목표 설정에 있어 어려움이 있다면 `SMART` 기법을 활용하는 것도 좋은 방법이다.

<img src="/md/img/unit-test-practice/smart-rule.png">
<em>S.M.A.R.T</em>

1. 테스트 케이스는 구체적이고 명확하게 작성한다.
2. 기능의 목표는 수치나 숫자로 측정할 수 있어야만 이후 문제가 있을 때 대비가 가능해진다.
3. 처음부터 무턱대고 기능의 난도를 높게 잡지 마라. 달성 가능한 목표를 세워야 한다.
4. 목표와 연관성이 있는 계획을 세워야 한다.
5. 모든 작업과 목표는 시간을 염두에 두고 수립되고 실행되어야 한다.

SMART 기법을 토대로 목표 설정을 정한 후 단위별로 나누어 테스트 케이스를 작성해보자.

<img src="/md/img/unit-test-practice/calcuration-test-case.png">
<em>문자열 계산기 테스트 케이스</em>

테스트 케이스를 작성했다. 테스트 케이스는 추후 개발을 하면서 수정 또는 추가할 수 있으므로 초반 테스트 케이스는 최대한 단순하고 명확하게 잡는 것이 중요하다. 

### 1. 테스트 클래스 생성한다. - Test.class Naming convention

 테스트 케이스에 첫 번째로 정의했던 테스트 클래스를 생성해보자.

_Calcuration.class -> CalcurationTest.class_
 
테스트 클래스의 이름은 테스트할 클래스 명 뒤에 Test를 붙여 생성한다. 이는 암묵적인 약속으로 테스트 클래스 명을 작성할 시 지켜야 할 암묵적인 명명규칙이다.

``` java
public class CalcurationTest {
    @Test
    public void test(){
    	Calcuration cal = new Calcuration();
    }
}
```

테스트 코드에 실제 코드가 작성될 `Calcuration` 클래스를 생성하였다. 다음 테스트 코드가 동작하는지 JUnit을 실행시켜보자.

<img src="/md/img/unit-test-practice/junit-execute.png">
<em>JUnit 실행 방법</em>

JUnit 실행 방법은 단축키 `Shift-Alt-X + T` 또는 컨텍스트 메뉴나 메인 메뉴에서 `Run As -> JUnit Test`로 실행할 수 있다.

<img src="/md/img/unit-test-practice/junit-gui-result1-1.png">
<em>[Red] Fail : Compile Error</em>

예상했듯이 `Calcuration cal = new Calcuration();` 부분이 컴파일 에러가 난다. 단위 테스트 작성에 있어 컴파일은 실패하지 않아야 한다는 규칙이 있다. 

_***규칙  1.Red 단계는 컴파일 에러가 아닌 코드상의 에러여야 한다.**_

컴파일 에러를 방지하기 위해 실제 구현할 src 소스 폴더에 Calcuration 클래스를 생성하자.

<img src="/md/img/unit-test-practice/class-structure.png">
<em>Calcuration 클래스 생성</em>

반드시 `Calcuration` 클래스 파일만 생성해야 하고 미리 앞서가서 메소드를 추가를 하면 안 된다. 이 말인즉슨 미리 앞서가서 다른 단계를 고려할 필요가 없다는 뜻이다.

_***규칙 2.테스트가 통과할 정도로만 실제 코드에 작성한다.**_

`Calcuration` 클래스를 생성했다면 다시 JUnit을 실행해보자. 

<img src="/md/img/unit-test-practice/junit-gui-result1-2.png">
<em>[Grean] Success</em>

이처럼 테스트 결과가 녹색이 나왔다면 테스트가 잘 진행되고 있다는 증거이다. 테스트 케이스에 체크 후 다음 단계를 진행하자.

> ~~1) 테스트 클래스를 생성한다.~~        <br/>
> 2) 문자열 입력값을 받는다.		  <br/>
> 3) 더하기 기능을 제공한다.        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

### 2. 문자열 입력값을 받는다.

문자열 입력값을 받아오는 방법은 다양하지만, `Calcuration` 생성자에 넘겨주고 get 메소드를 통해 입력값을 받아오는 방식을 선택했다. 물론 컴파일 에러는 나지 않도록 생성자와 get 메소드를 생성해주자.

``` java
public class CalcurationTest {
	@Test
	public void test(){
		String a = "1";
		String b = "2";
		
		Calcuration cal = new Calcuration(a, b);
		
		//assertEquals : 값을 비교하는 JUnit 단정문 메소드
		assertEquals(a, cal.getA());
		assertEquals(b, cal.getB());
	}
}

public class Calcuration {
	public Calcuration(String a, String b) {
	}

	public Object getA() {
		return null;
	}

	public Object getB() {
		return null;
	}
}
```

테스트 결과가 실패(Red)가 나오겠지만 JUnit을 실행해보자.

<img src="/md/img/unit-test-practice/junit-gui-result2-1.png">
<em>[Red] fail : java.lang.AssertionError</em>

예상했듯이 테스트 결과는 실패다. 어찌 보면 실패를 예상했으면서 코드를 추가하지 않고 실패를 기대하는 건 바보스러운 행동처럼 보인다. 

 하지만 이러한 행위는 TDD가 테스트 코드부터 작성하는 행위를 기반으로 두기 때문이다. 실제 코드를 먼저 작성하지 않고 테스트 코드를 먼저 작성하기 때문에 실패하는 단계를 보는 건 당연하다.
 
 이처럼 `Red` 단계를 기대하며 **실패하는 코드**를 <U>먼저 짜는 것은 매우 중요한 습관이다.</U> 코드를 검증하는 과정이 필요하므로, 꼭 Red 메시지를 활용하여 Green 단계에서 코드를 작성하지 않아 Red가 나왔다는 걸 알아야 한다. 
 
_***규칙 3.실제 코드를 작성하기 전 실패하는 테스트 코드부터 작성한다.**_

실패를 보았다면 앞서 설명한 `규칙 2`를 지키며  테스트 케이스를 통과할 정도로만 실제 코드를 수정하여  테스트를 통과해보자. 

``` java
public class Calcuration {
	int a;
	int b;
	
	public Calcuration(String a, String b) {
		this.a = Integer.parseInt(a);
		this.b = Integer.parseInt(b);
	}

	public int getA() {
		return a;
	}

	public int getB() {
		return b;
	}
}
```

녹색 바가 나왔다. 다음 단계를 진행하자.

<img src="/md/img/unit-test-practice/junit-gui-result2-2.png">
<em>[Grean] Success</em>

> ~~1) 테스트 클래스를 생성한다.~~        <br/>
> ~~2) 문자열 입력값을 받는다.~~  <br/>
> 3) 더하기 기능을 제공한다.        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

### 3. 더하기 기능을 제공한다.

cal.add() 메소드에 더하기 기능을 제공할 예정이다. `규칙1`을 준수하여 컴파일 에러를 방지하기 위해 실제 코드에도 add 메소드를 추가한다.  

``` java
public class CalcurationTest {
	...
	
	@Test
	public void addTest(){
		String a = "1";
		String b = "2";
		
		Calcuration cal = new Calcuration(a, b);
		assertEquals("더하기 테스트", 1+2, cal.add());
	}
}

public class Calcuration {
	...

	public int add(){
		return 0;
	}
}
```

JUnit을 실행하자. 당연히 실패하는 코드를 먼저 작성했기 때문에 테스트 결과는 실패다. 

<img src="/md/img/unit-test-practice/junit-gui-result3-1.png">
<em>[Red] fail : 더하기 테스트 ...</em>

먼저 실패하는 코드를 작성했다. 이후 테스트를 통과하기 위해 실제 코드를 수정 후 다시 테스트를 진행해보자. 

테스트가 성공적으로 진행됐다.

<img src="/md/img/unit-test-practice/junit-gui-result3-2.png">
<em>[Grean] Success</em>

#### RGR - 리펙토링(Clean Code)

앞서 단위 테스트 코드를 개발하는 데 있어 공통적인 주기가 나타난다.

이 주기는 테스트가 실패하는 `Red 단계`와 테스트가 성공하는 `Green 단계`를 반복하며 실제 코드를 작성하는 일련의 공통적인 주기가 나타났는데 이러한 반복적인 주기를 TDD-Cycle의 RGR 주기라고 한다. 

_RGR 주기 : Red[실패] -> Green[성공] -> Refactoring[리펙토링]_

특히나 실패와 성공 단계를 거쳤다면, 리펙토링 단계를 거쳐야 한다. 이 리펙토링 단계는 중복된 코드를 제거하여 최종적으로 `Clean Code`를 지향하기 때문이다. 

현재 우리는 더하기 기능만 추가했는데 벌써 중복된 코드가 보이기 시작했다. 테스트 코드를 리펙토링 하자.

``` java
public class CalcurationTest {

	String a = "1";
	String b = "2";
	
	Calcuration cal = new Calcuration(a, b);
	
	@Test
	public void test(){
		assertEquals("get a", a, cal.getA());
		assertEquals("get b", b, cal.getB());
	}
	
	@Test
	public void addTest(){
		assertEquals("더하기 테스트", 1+2, cal.add());
	}
}
```

 변수와 클래스 선언과 같은 중복된 코드를 최상단에 정의하여 리펙토링하였다. 
 
 반드시 리펙토링한 코드는 테스트 검증을 해야 한다. 테스트 클래스를 검증해보자.

<img src="/md/img/unit-test-practice/junit-gui-result4-1.png">
<em>[Grean] Success</em>

 성공적으로 테스트가 진행됐다. 리펙토링 단계에 의해 자연스레 코드 베이스가 깔끔해지고 명확해졌다. 이 과정은 `TDD`의 `ClaenCode` 장점을 간접적이나마 느낄 수 있는 과정이다.

**단 리팩토링한 코드는 반드시 통과되어야 한다.**

 테스트가 통과된다는 것은 로직상 이상이 없다는 것을 뜻하기 때문이다. 반면 리팩토링을 한 테스트 코드가 실패했다는 의미는 리팩토링이 잘못되었다는 의미와 직결한다.

마지막으로 실제 클래스도 리펙토링해보자.

``` java
public class Calcuration {
	int a;
	int b;
	
	public Calcuration(String a, String b) {
		this.a = Integer.parseInt(a);
		this.b = Integer.parseInt(b);
	}

	public int add() { return a + b; }
	
	public int getA() {
		return a;
	}

	public int getB() {
		return b;
	}
}
```

실제 클래스를 리펙토링 했다면 테스트 코드를 실행해보자. 테스트 결과가 성공적이다. 이처럼 테스트 클래스 뿐만 아니라 실제 클래스 또한 리펙토링을 하여 코드를 CleanCode 시켜줘야 한다.

### 4 ~ 6단계도 마찬가지로 같은 과정을 지키면서 실제코드에 작성한다.

> ~~1) 테스트 클래스를 생성한다.~~        <br/>
> ~~2) 문자열 입력값을 받는다.~~  	<br/>
> ~~3) 더하기 기능을 제공한다.~~      <br/>
> ~~4) 빼기 기능을 제공한다.~~        <br/>
> ~~5) 나누기 기능을 제공한다.~~       <br/>
> ~~6) 곱셈 기능을 제공한다.~~         <br/>

6단계를 끝으로 모든 테스트 케이스를 적용하여 단위 테스트 코드를 작성했다. 하지만 아직 문자열 계산기는 부족한 부분이 많다. 

>7) 입력 값이 null일때 처리

단편적인 예로 `입력 값이 null일때 처리` 같은 보완해야 할 로직이 많다. 이 경우엔 `TestCase를 작성하자`에서 설명했던 것처럼 테스트 케이스를 추가하여 점진적으로 기능을 보완하면 된다.

### 마무리

단위 테스트를 입문한 사람이라면 처음 코딩하는 습관이 매우 중요하다. 무엇보다 기능을 수정하거나 추가할 때마다 기능 검증을 통해 Red 메시지를 확인하는 습관이 매우 중요하다. 또한, 앞서 설명한 세 가지 단위 테스트 작성 규칙과 RGR 주기의 프로세스를 지키며 개발하여 습관을 기르는 것을 추천한다.

1. Red 단계는 컴파일 에러가 아닌 코드상의 에러여야 한다.
2. 테스트가 통과할 정도로만 실제 코드에 작성한다.
3. 실제 코드를 작성하기 전 실패하는 테스트 코드부터 작성한다.

<img src="/md/img/unit-test-practice/unit-test-flowchart.png">
<em>Unit Test 작성 순서도</em>

단위 테스트 작성은 기능의 설계에 대해 테스트 케이스를 작성하는 것부터 시작이다. 이 말인즉슨 설계가 어느 정도는 설계의 틀이 잡혀 있어야 하고 개발자가 전반적으로 기능에 대해 이해를 하고 있어야 한다는 뜻이다.

마지막으로 박재성님이 추천하는 연습 방법을 소개하며 글을 마무리한다.

#### 1단계 - xUnit 사용법과 단위 테스트 연습

1. 유틸성 메소드 테스트
2. 알고리즘 문제 테스트 
3. JDK 혹은 API의 사용법을 익히기 위한 학습 테스트

#### 2단계 - TDD 연습

1. 의존관계 (모바일 UI, 웹 UI, 데이터베이스)가 없고, 요구 사항이 명확한 프로그램으로 예제로 연습하기
2. 로또, 사다리타기, 볼링 점수판, 체스 등 게임류

#### 3단계 - 리팩토링과 객체지향 설계 또는 함수형 프로그래밍

1. 같은 요구사항을 반복해서 연습하기
2. 같은 과제를 구현할 때마다 제약 사항을 추가하여 난도를 높이면서 구현해야 실력이 향상된다.

>- 객체지향 설계 1단계 제약 사항
 1. 한 메소드에 오직 한단계의 들여쓰기만 한다.
 2. else를 쓰지 않는다.
- 객체지향 설계 1단계 제약 사항
 1. 모든 기본형 타입은 감싼다.
 2. 일급 컬렉션을 쓴다.
- 리팩토링 제약 사항
 1. 리팩토링 과정에서 컴파일 에러가 발생하지 않도록 한다.

#### 4단계 - 웹/모바일 프로젝트에 적용

1. 앞의 3단계와 병행한다.
2. 처음부터 바로 현장에 적용하는건 안된다.
3. 레거시에 적용하는건 진짜 난이도가 높은 기술이다.

#### 5단계 - ATDD와 CI 적용

1. ATDD (End to End 테스트)를 적용한다.
2. CI를 통해 지속적 피드백이 가능한 환경 구축한다.

### 참고

[UNIT Testing Tutorial - Learn in 10 Minutes](https://www.guru99.com/unit-testing-guide.html) <br/>
[What is a Unit Test?](https://content.pivotal.io/blog/what-is-a-unit-test-the-answer-might-surprise-you)<br/>
[백준 알고리즘](https://www.acmicpc.net/)<br/>
