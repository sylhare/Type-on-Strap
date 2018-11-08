---
layout: post
title: "TDD : TDD-Cycle 연습하기"
tags: [TDD, JUnit, JUnit4]
categories: [Test]
subtitle: "계산기 예제"
excerpt_separator: <!--more-->
display: "false"
sitemap:
changefreq: daily
priority: 1.0
---

<!--more-->

# TDD-Cycle 연습하기

---


### Test.class Naming convention

 가장 먼저 테스트 클래스를 생성하자.

_JUnitPractice.class -> JUnitPracticeTest.class_
 
테스트 클래스의 이름은 테스트할 클래스 명 뒤에 Test를 붙여 생성한다. 필수적이진 않지만, 일반적으로 자바 클래스의 명명규칙으로 첫 글자는 대문자로 짓는 것처럼 테스트 클래스 명을 지을 때 뒤에 Test를 붙여 작성한다. 이는 암묵적인 약속으로 테스트 클래스 명을 작성할 시 지켜야 할 암묵적인 명명규칙이다.

#### @Test public void methodName()

먼저 @Test 테스트 애노테이션은 메소드를 테스트 메소드로 정의한다. 앞서 '테스트를 수행하는 방식'에서 설명했듯이 JUnit이 테스트 메소드를 인식하기 위해선 3가지 규칙을 준수하여 메소드를 작성해야 한다.

1. 테스트할 메소드에 @Test 애노테이션을 지정한다.
2. 접근제한자는 반드시 public으로 정의한다.
3. 리턴 타입은 void로  테스트 메소드는 반환 값이 없다.

``` java
@Test public void test1(){ ... }
```

이처럼 간단한 작성 규칙만 따른다면 테스트를 할 수 있는 준비는 끝났다.

다음은 테스트 메소드의 결과를 검증해야될 차례다. 물론 검증에 있어  Self-Validating는 핵심이다. JUnit에서 지원하는 Assert 메소드를 통해 테스트를 자동화 검증을 할 수 있다.

>단위 테스트의 FIRST 규칙 중 Self-Validating 규칙은 수동 테스트가 아닌 자동화 테스트가 되어야 한다는 의미로 테스트 결과가 올바른지에 대한 판단은 개발자가 임의로 결정해서는 안 된다. 

### TDD 연습하기

TDD를 개발 연습할 때는 유틸성 기능 또는 알고리즘 문제를 통해 연습하면 좋다.<br/>
유틸성 문제 중 가장 기초적인 문자열 계산기를 만들어보자.<br/>

> 문제 1) 더하기, 빼기, 나누기, 곱셈할 수 있는 문자열 계산기 만들기


#### 먼저 TestCase를 작성하자

`TDD`는 테스트 코드를 작성하기 이전에 테스트 케이스를 추출하는 것부터 시작한다.<br/>
개발자는 테스트 케이스를 통해 개발의 진척도를 직관적으로 알 수 있다.
또한, 목록을 하나씩 지울 때마다 소소한 성취감을 느낄 수 있고 무엇보다 해야 할 일들이 눈에 보이기 때문에 앞만 보고 개발을 진행할 수 있다.

##### 목표 세우기

"테스트 케이스를 어떻게 작성하면 좋을까?"라는 의문이 생긴다. 이러한 문제는 목표설정에 있어 `SMART` 기법을 활용하여 편하게 접근할 수 있다.<br/>

<img src="/md/img/tdd-practice/smart-rule.png" height="400px">
<em>S.M.A.R.T</em>

이후 목표를 정한 후 단위별로 나누어 테스트 케이스를 작성한다.

> **계산기 테스트 케이스** 			<br/>
> 1) 계산기 클래스 생성           		<br/>
> 2) 문자열 입력값을 받는다.  <br/>
> 3) 더하기 기능을 제공한다.        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

1) 계산기 테스트 클래스를 생성한다. <br/>
   - 테스트 클래스의 클래스 명은 실제 구현할 클래스 뒤에 Test를 붙여 만든다. (암묵적인 규칙이다.)


```java
import org.junit.Test;

public class CalcurationTest {

    @Test
    public void test(){
    	Calcuration cal = new Calcuration();
    }
}
```

*'컴파일은 실패하지 않으면서 실행이 실패하는 정도로만 테스트 코드를 작성한다. '*

위의 코드는 실제 코드가 작성될 `Calcuration` 클래스를 생성하였다. <br/>
당연하지만 `Calcuration cal = new Calcuration();` 부분이 에러가 난다. 오류 해결을 위해 `Calcuration` 클래스 파일을 생성한다. <br/>
반드시 `Calcuration` 클래스 파일만 생성하자.(미리 앞서가서 메소드를 추가를 하면 안 된다.)<br/>
`Calcuration` 클래스를 생성한 후 테스트를 해보자. 테스트가 녹색 바가 나왔다면 테스트가 잘 진행되고 있다는 증거이다. 테스트 케이스에 체크 후 다음 단계를 진행하자.<br/>

> **계산기 테스트 케이스** 			<br/>
> ~~1) 계산기 클래스 생성~~           		<br/>
> 2) 문자열 입력값을 받는다.  <br/>
> 3) 더하기 기능을 제공한다.        <br/>
> 4) 빼기 기능을 제공한다.         <br/>
> 5) 나누기 기능을 제공한다.        <br/>
> 6) 곱셈 기능을 제공한다.         <br/>

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

---
