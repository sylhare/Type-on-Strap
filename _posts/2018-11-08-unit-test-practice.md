---
layout: post
title: "단위 테스트 코드를 작성 해보자"
tags: [TDD, TDDCycle, RGRCycle, JUnit, JUnit4, UnitTest]
categories: [Test]
subtitle: "RGR 단계를 인식하며 단위 테스트 코드를 개발하기"
feature-img: "md/img/thumbnail/tdd-cycle-practice.png"    
thumbnail: "md/img/thumbnail/tdd-cycle-practice.png"
excerpt_separator: <!--more-->
display: "false"
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

이처럼 정의한 테스트 케이스는 추후 개발을 하면서 수정할 수 있으므로 초반 테스트 케이스는 최대한 단순하고 명확하게 잡는 것이 중요하다. 

### 1. 테스트 클래스 생성한다. - Test.class Naming convention

 테스트 케이스에 첫 번째로 정의했던 테스트 클래스를 생성해보자.

_Calcuration.class -> CalcurationTest.class_
 
테스트 클래스의 이름은 테스트할 클래스 명 뒤에 Test를 붙여 생성한다. 이는 암묵적인 약속으로 테스트 클래스 명을 작성할 시 지켜야 할 암묵적인 명명규칙이다.

```java
import org.junit.Test;

public class CalcurationTest {

    @Test
    public void test(){
    	Calcuration cal = new Calcuration();
    }
}
```

테스트 코드에 실제 코드가 작성될 `Calcuration` 클래스를 생성하였다. 다음 테스트 코드가 동작하는지 JUnit을 실행시켜보자.

<img src="/md/img/unit-test-practice/junit-execute.png">
<em>Run As -> JUnit Test</em>

JUnit 실행은 단축키 `Shift-Alt-X + T` 또는 컨텍스트 메뉴나 메인 메뉴에서 `Run As -> JUnit Test`로 실행할 수 있다.

<img src="/md/img/unit-test-practice/junit-gui-result1-1.png">
<em>[Red] Fail : Compile Error</em>

예상했듯이 `Calcuration cal = new Calcuration();` 부분이 컴파일 에러가 난다. 단위 테스트 작성에 있어 컴파일은 실패하지 않아야 한다는 규칙이 있다. 

*규칙 1. 컴파일은 실패하지 않으면서 실행이 실패하는 정도로만 테스트 코드를 작성한다.*

컴파일 에러를 방지하기 위해 실제 구현할 src 소스 폴더에 Calcuration 클래스 또한 생성해주자.

<img src="/md/img/unit-test-practice/class-structure.png">
<em>Calcuration 클래스 생성</em>

반드시 `Calcuration` 클래스 파일만 생성하자. 미리 앞서가서 메소드를 추가를 하면 안 된다. 단위 테스트는 한 단계를 지키며 나아가는게 중요하기 때문이다.

`Calcuration` 클래스를 생성했다면 다시 JUnit을 실행해보자. 

<img src="/md/img/unit-test-practice/junit-gui-result1-2.png">
<em>[Grean] Success</em>

테스트가 녹색 바가 나왔다면 테스트가 잘 진행되고 있다는 증거이다. 테스트 케이스에 체크 후 다음 단계를 진행하자.<br/>

~~1. 테스트 클래스 생성~~ 

---

3) 두 개의 정수형인 입력값 받는다.

##### CalcurationTest.java
```java
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CalcurationTest {
	
	@Test
	public void test(){
		String a = "1";
		String b = "2";
		
		Calcuration cal = new Calcuration(a, b);
		
		assertEquals("get a", "1", cal.getA(), 0);
		assertEquals("get b", "2", cal.getB(), 0);
	}
	
}
```

`Calcuration` 생성자에 정수 a와 b를 넘겨주고 get 메소드를 통해 입력 값이 잘 전달됐는지 확인하는 테스트 코드이다. <br/>
정숫값 비교를 위해 `assertEquals` 메소드를 사용했고 사용법은  `assertEquals( 기대하는 값, 값, 오차범위);`이다.<br/>
테스트하면 당연히 실패이다. 실제 코드에 생성자와 get 메소드를 추가하자.

##### Calcuration.java

``` java
package com.test.jUnitPratice.tddJunit;

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

*'TDD는 반드시 현재 실패하는 테스트 케이스를 통과할 정도로만 실제 코드에 작성한다.'*
실제 코드에 작성한 후 테스트를 진행해보자.
녹색 바가 나왔다. 다음 단계를 진행하자.

> **계산기 테스트 케이스** 			<br/>
> ~~1) 계산기 클래스 생성~~        <br/>
> ~~2) 문자열 입력값을 받는다.~~  <br/>
> 3) 더하기 기능을 제공한다.        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

3) 더하기 기능을 제공한다.

##### CalcurationTest.java
``` java
@Test
public void add(){
	String a = "1";
	String b = "2";
	
	Calcuration cal = new Calcuration(a, b);
	
	assertEquals("더하기 테스트", 1+2, cal.add());
}
```

실패하는 코드를 작성한다. 테스트를 통과하기 위해 실제 코드에 add 메소드를 작성 후 다시 테스트를 진행해보자.<br/>
테스트가 성공적으로 진행됐다.

#### 자연스럽게 보이는 RGR 주기 (ClaenCode)

*'실패 -> 성공 -> 리펙토링  (RGR 주기)'*

더하기 기능만 추가했는데 벌써 중복된 코드가 보이기 시작한다. <br/>
테스트 코드를 리펙토링 후 테스트를 진행해보자.

##### CalcurationTest.java
``` java
package com.test.jUnitPratice.tddJunit;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class CalcurationTest {

	String a = "1";
	String b = "2";
	
	Calcuration cal = new Calcuration(a, b);
	
	@Test
	public void test(){
		assertEquals("get a", 1, cal.getA(), 0);
		assertEquals("get b", 2, cal.getB(), 0);
	}
	
	@Test
	public void add(){
		assertEquals("더하기 테스트", 1+2, cal.add());
	}
}
```

성공적으로 테스트가 진행됐다. `RGR` 주기에 의해 자연스레 코드 베이스가 깔끔해지고 명확해진다.
이는 `TDD`의 `ClaenCode` 장점을 간접적이나마 느낄 수 있는 과정이다. <br/>

단 리팩토링한 코드는 반드시 통과되어야 한다. 테스트가 통과된다는 것은 로직상 이상이 없다는 것을 뜻하기 때문이다.<br/>
리팩토링을 한 테스트 코드가 실패했다는 의미는 리팩토링이 잘못되었다는 의미와 직결한다.



> **계산기 테스트 케이스** 			<br/>
> ~~1) 계산기 클래스 생성~~        <br/>
> ~~2) 문자열 입력값을 받는다.~~  <br/>
> ~~3) 더하기 기능을 제공한다.~~        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

---

### 4 ~ 6단계도 마찬가지로 같은 과정을 지키면서 실제코드에 작성한다.

##### Calcuration.java
``` java
public class Calcuration {
	
	int a;
	int b;
	
	public Calcuration(String a, String b) {
		this.a = Integer.parseInt(a);
		this.b = Integer.parseInt(b);
	}


	public int add() { return a + b; }
	public int sub() { return a - b; }
	public int muliply() { return a * b; }
	public int division() { return a % b; }

	public int getA() {
		return a;
	}

	public int getB() {
		return b;
	}
}
```

##### CalcurationTest.java
``` java
import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class CalcurationTest {
	
	String a = "1";
	String b = "2";
	
	Calcuration cal = new Calcuration(a, b);
	
	@Test
	public void test(){
		assertEquals("get a", Integer.parseInt(a), cal.getA(), 0);
		assertEquals("get b", Integer.parseInt(b), cal.getB(), 0);
	}
	
	@Test
	public void add(){
		assertEquals("더하기 테스트", Integer.parseInt(a)+Integer.parseInt(b), cal.add());
	}
	
	@Test
	public void sub(){
		assertEquals("빼기 테스트", Integer.parseInt(a)-Integer.parseInt(b), cal.sub());
	}
	
	@Test
	public void muliply(){
		assertEquals("곱하기 테스트", Integer.parseInt(a)*Integer.parseInt(b), cal.muliply());
	}
	
	@Test
	public void division(){
		assertEquals("나누기 테스트", Integer.parseInt(a)%Integer.parseInt(b), cal.division());
	}
}
```

---

### Goal?

> **계산기 테스트 케이스** 			<br/>
> ~~1) 계산기 클래스 생성~~        <br/>
> ~~2) 문자열 입력값을 받는다.~~  <br/>
> ~~3) 더하기 기능을 제공한다.~~        <br/>
> ~~4) 빼기 기능을 제공한다.~~         <br/>
> ~~5) 나누기 기능을 제공한다.~~        <br/>
> ~~6) 곱셈 기능을 제공한다.~~         <br/>

테스트 케이스를 끝냈다. 

### 마무리

여기까지 `블랙 박스 테스트`와 `화이트 박스 테스트`를 통해 로직을 검증하는 과정에서 발견되는 여러 에러 사항을 처리함으로써 로직을 보완했다.
`블랙 박스 테스트`에서는 입출력에 대한 검증을 통해 초기 데이터의 에러에 대해 결함률을 낮추고 `화이트 박스 테스트`에서는  오류 주입을 통해 예외 처리 코드의 분기에 사용되어 로직을 보완했다.<br/>

여러 단위 테스트를 할 수 있는 `Library` 중에 `JUnit`을 선택한 이유는 <br/>
첫 번째 나는 자바 개발자이기 때문이다. <br/>
두 번째 `JUnit`은 `Java`의 테스트 `Library` 중 선두주자이기 때문이다.(5버전이 나온 이유다.) <br/>
세 번째 `Spring boot`환경에서 `JUnit5`을 접목해 테스트하는 방법이 최근에 나오기 시작했다. <br/>
마지막으론 기술 트렌드를 놓치고 싶지 않기 위함이다. 추후 `Spring boot`에서 `JUnit5`을 접목해 테스트하는 방법을 작성할 계획이다.

### 참고

[https://www.guru99.com/unit-testing-guide.html](https://www.guru99.com/unit-testing-guide.html) <br/>
[https://content.pivotal.io/blog/what-is-a-unit-test-the-answer-might-surprise-you]<br/>


https://jojoldu.tistory.com/306




---
